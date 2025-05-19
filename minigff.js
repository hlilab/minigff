#!/usr/bin/env k8

"use strict";

const gff_version = "r12";

/*********************************
 * Command-line argument parsing *
 *********************************/

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

/******************
 * Interval query *
 ******************/

function iit_sort_dedup_copy(a) {
	a.sort((x, y) => (x.st != y.st? x.st - y.st : x.en - y.en));
	const b = [];
	for (let i0 = 0, i = 1; i <= a.length; ++i) {
		if (i == a.length || a[i0].st != a[i].st || a[i0].en != a[i].en) {
			b.push({ st:a[i0].st, en:a[i0].en, max:0 });
			i0 = i;
		}
	}
	return b;
}

function iit_index(a) {
	if (a.length == 0) return -1;
	let last, last_i, k;
	for (let i = 0; i < a.length; i += 2) last = a[i].max = a[i].en, last_i = i;
	for (k = 1; 1<<k <= a.length; ++k) {
		const i0 = (1<<k) - 1, step = 1<<(k+1), x = 1<<(k-1);
		for (let i = i0; i < a.length; i += step) {
			a[i].max = a[i].en;
			if (a[i].max < a[i-x].max) a[i].max = a[i-x].max;
			const e = i + x < a.length? a[i+x].max : last;
			if (a[i].max < e) a[i].max = e;
		}
		last_i = last_i>>k&1? last_i - x : last_i + x;
		if (last_i < a.length) last = last > a[last_i].max? last : a[last_i].max;
	}
	return k - 1;
}

function iit_overlap(a, st, en) {
	let h = 0;
	const stack = [], b = [];
	for (h = 0; 1<<h <= a.length; ++h);
	--h;
	stack.push([(1<<h) - 1, h, 0]);
	while (stack.length) {
		const t = stack.pop();
		const x = t[0], h = t[1], w = t[2];
		if (h <= 3) {
			const i0 = x >> h << h;
			let i1 = i0 + (1<<(h+1)) - 1;
			if (i1 >= a.length) i1 = a.length;
			for (let i = i0; i < i1 && a[i].st < en; ++i)
				if (st < a[i].en) b.push(a[i]);
		} else if (w == 0) { // if left child not processed
			stack.push([x, h, 1]);
			const y = x - (1<<(h-1));
			if (y >= a.length || a[y].max > st)
				stack.push([y, h - 1, 0]);
		} else if (x < a.length && a[x].st < en) {
			if (st < a[x].en) b.push(a[x]);
			stack.push([x + (1<<(h-1)), h - 1, 0]);
		}
	}
	return b;
}

/********************
 * Simpler File I/O *
 ********************/

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

/***************
 *** minigff ***
 ***************/

/******************
 * Format parsing *
 ******************/

class Transcript {
	constructor(tid, ctg, strand) {
		this.tid = tid, this.ctg = ctg, this.strand = strand;
		this.st = this.en = this.cds_st = this.cds_en = -1;
		this.disp_name = null;
		this.type = null;
		this.gid = null, this.gname = null;
		this.pri = null;
		this.score = -1;
		this.err = 0, this.done = false;
		this.exon = [];
	}
	finish() {
		if (this.exon.length == 0) return false;
		const a = this.exon.sort(function(a,b) {return a.st - b.st});
		this.exon = a;
		const st = a[0].st;
		const en = a[a.length - 1].en;
		if (this.st < 0) this.st = st;
		if (this.en < 0) this.en = en;
		if (st != this.st || en != this.en) this.err |= 1;
		for (let i = 1; i < a.length; ++i)
			if (a[i].st < a[i-1].en)
				this.err |= 2;
		if ((this.cds_st >= 0 && this.cds_st < st) || (this.cds_en >= 0 && this.cds_en > en)) this.err |= 4;
		if (!this.disp_name)
			this.disp_name = this.gname && this.type? [this.tid, this.type, this.gname].join("|") : this.tid;
		if (this.err == 0) this.done = true;
		if (!this.done) this.perror();
		return this.done;
	}
	perror() {
		if (this.err & 1) warn(`ERROR: inconsistent transcript and exon coordinates for transcript ${this.tid}`);
		if (this.err & 2) warn(`ERROR: overlapping exon for transcript ${this.tid}`);
		if (this.err & 4) warn(`ERROR: inconsistent thick start or end for transcript ${this.tid}`);
	}
	bed() {
		let offs = [], lens = [];
		for (let i = 0; i < this.exon.length; ++i) {
			lens.push(this.exon[i].en - this.exon[i].st);
			offs.push(this.exon[i].st - this.st);
		}
		return [this.ctg, this.st, this.en, this.disp_name, this.score < 0? "." : this.score, this.strand,
			this.cds_st < 0? "." : this.cds_st, this.cds_en < 0? "." : this.cds_en, ".",
			this.exon.length, lens.join(",")+",", offs.join(",")+","];
	}
}

function* gff_read(fn, cds_only, disp_target_name) {
	const re_fmt = /^(##PAF\t\S+(\t\d+){3}\t[+-])|^(\S+(\t\d+){3}\t[+-])|^(\S+\t\d+\t\d+\t\S+\t\S+\t[+-])|^((\S+\t){3}\d+\t\d+\t\S+\t[+-])|^(\S+\t\d+\t\S+\t\d+\t\d+\t(\d+[MIDNSH=X])+)/;
	const re_gff = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|tag|Parent|Target)( "([^"]+)"|=([^;]+))/g;
	const re_cigar = /(\d+)([MIDNSHP=XFGUV])/g;
	if (typeof cds_only == "undefined") cds_only = false;
	if (typeof disp_target_name == "undefined") disp_target_name = false;

	function nt_cigar2exon(v, st, cigar) {
		let m, x = st, x0 = st;
		while ((m = re_cigar.exec(cigar)) != null) {
			const len = parseInt(m[1]), op = m[2];
			if (op == 'N') {
				v.exon.push({ st:x0, en:x });
				x0 = x + len, x += len;
			} else if (op == 'M' || op == 'X' || op == '=' || op == 'D') {
				x += len;
			}
		}
		v.exon.push({ st:x0, en:x });
		return x;
	}

	let v = null, fmt0 = -1;
	for (const line of k8_readline(fn)) {
		let m;
		if ((m = re_fmt.exec(line)) == null) continue;
		let fmt = m[5]? 1 : m[6]? 2 : m[3]? 3 : m[1]? 4 : m[8]? 5 : -1; // 1=BED, 2=GFF, 3=PAF, 4=##PAF, 5=SAM
		if (fmt < 0) continue;
		if (fmt0 < 0) fmt0 = fmt;
		if (fmt0 != fmt) continue; // only allow one format
		if (fmt == 1) { // BED12
			let t = line.split("\t");
			if (t.length < 12) continue;
			if (cds_only && (t[6] == "." || t[7] == ".")) continue;
			v = new Transcript(t[3], t[0], t[5]);
			v.st = parseInt(t[1]);
			v.en = parseInt(t[2]);
			v.cds_st = t[6] == "."? -1 : parseInt(t[6]);
			v.cds_en = t[7] == "."? -1 : parseInt(t[7]);
			const n_exon = parseInt(t[9]);
			const lens = t[10].split(",", n_exon);
			const offs = t[11].split(",", n_exon);
			for (let i = 0; i < n_exon; ++i) {
				const off = parseInt(offs[i]), len = parseInt(lens[i]);
				v.exon.push({ st: v.st + off, en: v.st + off + len });
			}
			if (v.finish()) yield v;
			v = null;
		} else if (fmt == 2) { // GTF or GFF3
			let t = line.split("\t");
			if (t[2] != "exon" && t[2] != "CDS" && t[2] != "cds") continue;
			let m, tid = null, gid = null, type = "", gname = null, biotype = "", ens_canonical = false, pid = null, target = null;
			while ((m = re_gff.exec(t[8])) != null) {
				const key = m[1], val = m[3]? m[3] : m[4];
				if (key == "transcript_id") tid = val;
				else if (key == "transcript_type") type = val;
				else if (key == "transcript_biotype" || key == "gbkey") biotype = val;
				else if (key == "gene_name") gname = val;
				else if (key == "gene_id") gid = val;
				else if (key == "Parent") pid = val;
				else if (key == "tag" && val == "Ensembl_canonical") ens_canonical = true;
				else if (key == "Target") target = val.split(" ")[0];
			}
			if (tid == null) tid = pid;
			if (tid == null) continue;
			if (gname == null) gname = gid? gid : "*"; // infer gene name
			if (gid == null) gid = gname; // if gene_id missing, use gene name to identify a gene
			if (type == "" && biotype != "") type = biotype; // infer transcript type
			if (v == null || v.tid != tid) {
				if (v != null && v.finish()) yield v;
				v = new Transcript(tid, t[0], t[6]);
				v.type = type, v.gid = gid, v.gname = gname;
			}
			if (disp_target_name && target != null)
				v.disp_name = target;
			const st = parseInt(t[3]) - 1;
			const en = parseInt(t[4]);
			if (t[2] === "CDS" || t[2] === "cds") {
				v.cds_st = v.cds_st >= 0 && v.cds_st < st? v.cds_st : st;
				v.cds_en = v.cds_en >= 0 && v.cds_en > en? v.cds_en : en;
				if (cds_only) v.exon.push({ st:st, en:en });
			} else if (t[2] === "exon" && !cds_only) {
				v.exon.push({ st:st, en:en });
			}
		} else if (fmt == 3 || fmt == 4) { // PAF
			let m, t = line.split("\t");
			if (fmt == 4) t.shift();
			let cigar = null, pri = null, ms = null, score = null;
			for (let i = 12; i < t.length; ++i) {
				const key = t[i].substr(0, 5);
				if (key == "cg:Z:") cigar = t[i].substr(5);
				else if (key == "tp:A:") pri = t[i].substr(5) == "P"? true : false;
				else if (key == "ms:i:") ms = parseInt(t[i].substr(5));
				else if (key == "AS:i:") score = parseInt(t[i].substr(5));
			}
			if (cigar == null) continue;
			if (ms != null) score = ms;
			const st = parseInt(t[7]);
			const en = parseInt(t[8]);
			let rlen = 0, is_aa = false;
			while ((m = re_cigar.exec(cigar)) != null) {
				const op = m[2], len = parseInt(m[1]);
				if (op == "U" || op == "V" || op == "F" || op == "G") is_aa = true;
				else if (op == "M" || op == "X" || op == "=" || op == "N" || op == "D") rlen += len;
			}
			if (!is_aa)
				is_aa = rlen == en - st? false : true;
			v = new Transcript(t[0], t[5], t[4]);
			if (score != null) v.score = score;
			if (pri != null) v.pri = pri;
			if (is_aa) { // amino acid PAF
				let x = 0, x0 = 0, e = [];
				while ((m = re_cigar.exec(cigar)) != null) {
					const len = parseInt(m[1]), op = m[2];
					if (op == 'N') {
						e.push([x0, x]);
						x0 = x + len, x += len;
					} else if (op == 'U') {
						e.push([x0, x + 1]);
						x0 = x + len - 2, x += len;
					} else if (op == 'V') {
						e.push([x0, x + 2]);
						x0 = x + len - 1, x += len;
					} else if (op == 'M' || op == 'X' || op == '=' || op == 'D') {
						x += len * 3;
					} else if (op == 'F' || op == 'G') {
						x += len;
					}
				}
				e.push([x0, x]);
				if (x != en - st) throw Error("inconsistent CIGAR");
				if (t[4] == '+') {
					for (let i = 0; i < e.length; ++i)
						v.exon.push({ st: st + e[i][0], en: st + e[i][1] });
				} else if (t[4] == '-') { // For protein-to-genome alignment, the coordinates are on the query strand. Need to flip them.
					const glen = en - st;
					for (let i = e.length - 1; i >= 0; --i)
						v.exon.push({ st: st + (glen - e[i][1]), en: st + (glen - e[i][0]) });
				}
			} else { // nucleotide PAF
				let x = nt_cigar2exon(v, st, cigar);
				if (x != en) throw Error("inconsistent CIGAR");
			}
			if (v.finish()) yield v;
			v = null;
		} else if (fmt == 5) { // SAM
			let m, t = line.split("\t");
			let ms = null, score = null;
			for (let i = 11; i < t.length; ++i) {
				const key = t[i].substr(0, 5);
				if (key == "ms:i:") ms = parseInt(t[i].substr(5));
				else if (key == "AS:i:") score = parseInt(t[i].substr(5));
			}
			if (ms != null) score = ms;
			const flag = parseInt(t[1]);
			v = new Transcript(t[0], t[2], flag&16? "-" : "+");
			v.score = score;
			v.pri = flag&0x100? false : true;
			nt_cigar2exon(v, parseInt(t[3]) - 1, t[5]);
			if (v.finish()) yield v;
			v = null;
		}
	}
	if (v != null && v.finish()) yield v;
}

function gff_cmd_all2bed(args)
{
	let pri_only = false, cds_only = false, disp_target_name = false, print_junc = false;
	for (const o of getopt(args, "aptj", [])) {
		if (o.opt == "-p") pri_only = true;
		else if (o.opt == "-a") cds_only = true;
		else if (o.opt == "-t") disp_target_name = true;
		else if (o.opt == "-j") print_junc = true;
	}
	if (args.length == 0) {
		print("Usage: minigff.js all2bed [options] <in.file>");
		print("Options:");
		print("  -a       only process CDS");
		print("  -p       only include primary alignments");
		print("  -j       print junctions/introns");
		print("  -t       display Target name");
		return;
	}
	for (let v of gff_read(args[0], cds_only, disp_target_name)) {
		if (pri_only && !v.pri) continue;
		if (print_junc) {
			for (let i = 1; i < v.exon.length; ++i)
				print(v.ctg, v.exon[i-1].en, v.exon[i].st, ".", ".", v.strand);
		} else {
			print(v.bed().join("\t"));
		}
	}
}

/************
 * Evaluate *
 ************/

class BaseIndex {
	constructor() {
		this.ctg = {};
	}
	add(v, cds_only) { // add a Transcript object
		if (this.ctg[v.ctg] == null) this.ctg[v.ctg] = { exon:[], junc:[] };
		let p = this.ctg[v.ctg], exon;
		if (cds_only) {
			exon = [];
			for (let i = 0; i < v.exon.length; ++i) {
				const st = v.cds_st > v.exon[i].st? v.cds_st : v.exon[i].st;
				const en = v.cds_en < v.exon[i].en? v.cds_en : v.exon[i].en;
				if (st < en) exon.push({ st:st, en:en });
			}
		} else exon = v.exon;
		for (let i = 0; i < exon.length; ++i)
			p.exon.push({ st:exon[i].st, en:exon[i].en });
		for (let i = 1; i < exon.length; ++i)
			p.junc.push({ st:exon[i-1].en, en:exon[i].st });
	}
	index() {
		for (const name in this.ctg) {
			let p = this.ctg[name];
			p.exon = iit_sort_dedup_copy(p.exon);
			p.junc = iit_sort_dedup_copy(p.junc);
			iit_index(p.exon);
			iit_index(p.junc);
		}
	}
	ov_exon(name, st, en) {
		return this.ctg[name] == null? [] : iit_overlap(this.ctg[name].exon, st, en);
	}
	ov_junc(name, st, en) {
		return this.ctg[name] == null? [] : iit_overlap(this.ctg[name].junc, st, en);
	}
}

function gff_cmd_eval(args)
{
	let print_all = false, print_err = false, first_only = false, chr_only = false, skip_last = false, skip_first = false, cds_only = false, eval_base = false;
	for (const o of getopt(args, "e1ctfdaps", [])) {
		if (o.opt == "-e") print_err = true;
		else if (o.opt == "-p") print_all = true;
		else if (o.opt == "-1") first_only = true;
		else if (o.opt == "-c") chr_only = true;
		else if (o.opt == "-f") skip_first = true;
		else if (o.opt == "-t") skip_first = skip_last = true;
		else if (o.opt == "-s") eval_base = true;
		else if (o.opt == "-d" || o.opt == "-a") cds_only = true;
	}
	if (args.length < 2) {
		print("Usage: minigff.js eval [options] <base.file> <test.file>");
		print("Options:");
		print("  -a      ignore UTRs in BASE");
		print("  -1      only evaluate the first alignment of each TEST");
		print("  -c      only consider TEST alignments to contig /^(chr)?([0-9]+|X|Y)$/");
		print("  -f      skip the first exon in TEST for exon evaluation");
		print("  -t      skip the last exon in TEST for exon evaluation");
		print("  -s      evaluate base Sn and Sp (more memory)");
		print("  -e      print error intervals");
		print("  -p      print all intervals");
		return;
	}

	function gff_format_ov(ov) {
		let x = '[';
		for (let j = 0; j < ov.length; ++j) {
			if (j) x += ', ';
			x += `(${ov[j].st},${ov[j].en})`
		}
		x += ']';
		return x;
	}

	print("CC\tNN  nTest  nSingleton");
	print("CC\tNE  nExon  nAnnoExon  ratio  nNovelExon");
	print("CC\tNJ  nJunc  nAnnoJunc  ratio  nNovelJunc");
	print("CC");

	// load base annotation
	let base = new BaseIndex();
	for (let v of gff_read(args[0], cds_only))
		base.add(v, cds_only);
	base.index();

	// evaluate test annotation
	let last_tid = null, tot_exon = 0, ann_exon = 0, nov_exon = 0, tot_junc = 0, ann_junc = 0, nov_junc = 0, n_multi = 0, n_test = 0;
	let qexon = {};
	for (let v of gff_read(args[1])) {
		if (chr_only && !/^(chr)?([0-9]+|X|Y)$/.test(v.ctg)) continue;
		if (first_only && last_tid == v.tid) continue;
		last_tid = v.tid;
		++n_test;
		if (v.exon.length > 1) ++n_multi;
		// skip first/last exon if requested
		let i_st = 0, i_en = v.exon.length;
		if (v.strand == "+") {
			if (skip_first) ++i_st;
			if (skip_last)  --i_en;
		} else if (v.strand == "-") {
			if (skip_first) --i_en;
			if (skip_last)  ++i_st;
		}
		for (let i = i_st; i < i_en; ++i) { // test exons
			const st = v.exon[i].st, en = v.exon[i].en;
			let found = 0, ov = base.ov_exon(v.ctg, st, en);
			if (ov.length == 0) ++nov_exon;
			for (let j = 0; j < ov.length; ++j)
				if (ov[j].st == st && ov[j].en == en)
					++found;
			++tot_exon;
			if (found > 0) ++ann_exon;
			if (eval_base) {
				if (qexon[v.ctg] == null) qexon[v.ctg] = [];
				qexon[v.ctg].push({ st:st, en:en });
			}
			if (print_all || (found == 0 && print_err)) {
				const label = ov.length == 0? "EN" : found == 0? "EP" : "EC";
				if (ov.length > 0) print(label, v.tid, i+1, v.ctg, st, en, gff_format_ov(ov));
				else print(label, v.tid, i+1, v.ctg, st, en);
			}
		}
		for (let i = 1; i < v.exon.length; ++i) { // test junctions
			const st = v.exon[i-1].en, en = v.exon[i].st;
			let found = 0, ov = base.ov_junc(v.ctg, st, en);
			if (ov.length == 0) ++nov_junc;
			for (let j = 0; j < ov.length; ++j)
				if (ov[j].st == st && ov[j].en == en)
					++found;
			++tot_junc;
			if (found > 0) ++ann_junc;
			if (print_all || (found == 0 && print_err)) {
				const label = ov.length == 0? "JN" : found == 0? "JP" : "JC";
				if (ov.length > 0) print(label, v.tid, i, v.ctg, st, en, gff_format_ov(ov));
				else print(label, v.tid, i, v.ctg, st, en);
			}
		}
	}
	print("NN", n_test, n_test - n_multi);
	print("NE", tot_exon, ann_exon, (ann_exon / tot_exon).toFixed(4), nov_exon);
	print("NJ", tot_junc, ann_junc, (ann_junc / tot_junc).toFixed(4), nov_junc);

	// compute base Sn/Sp
	function merge_and_index(ex) {
		for (const chr in ex) {
			let a = [], e = iit_sort_dedup_copy(ex[chr]);
			let st = e[0].st, en = e[0].en;
			for (let i = 1; i < e.length; ++i) { // merge
				if (e[i].st > en) {
					a.push({ st:st, en:en });
					st = e[i].st, en = e[i].en;
				} else {
					en = en > e[i].en? en : e[i].en;
				}
			}
			a.push({ st:st, en:en });
			iit_index(a);
			ex[chr] = a;
		}
	}

	function cal_sn(a0, a1) {
		let tot = 0, cov = 0;
		for (const chr in a1) {
			let e0 = a0[chr], e1 = a1[chr];
			for (let i = 0; i < e1.length; ++i)
				tot += e1[i].en - e1[i].st;
			if (e0 == null) continue;
			for (let i = 0; i < e1.length; ++i) {
				const o = iit_overlap(e0, e1[i].st, e1[i].en);
				for (let j = 0; j < o.length; ++j) { // this only works when there are no overlaps between intervals
					let st = e1[i].st > o[j].st? e1[i].st : o[j].st;
					let en = e1[i].en < o[j].en? e1[i].en : o[j].en;
					cov += en - st;
				}
			}
		}
		return [tot, cov];
	}

	if (eval_base) {
		let bexon = {};
		for (const ctg in base.ctg)
			bexon[ctg] = base.ctg[ctg].exon;
		merge_and_index(qexon);
		merge_and_index(bexon);
		const sn = cal_sn(qexon, bexon);
		const sp = cal_sn(bexon, qexon);
		print("BN", sn[0], sn[1], (sn[1] / sn[0]).toFixed(4));
		print("BP", sp[0], sp[1], (sp[1] / sp[0]).toFixed(4));
	}
}

/*****************
 * Main function *
 *****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: minigff.js <command> [arguments]");
		print("Commands:");
		print("  all2bed        convert BED12/GFF/GTF/PAF/SAM to BED12");
		print("  eval           evaluate against reference annotations");
		print("  version        print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == "all2bed") gff_cmd_all2bed(args);
	else if (cmd == "eval") gff_cmd_eval(args);
	else if (cmd == "version") {
		print(gff_version);
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
