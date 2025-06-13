#!/usr/bin/env k8

"use strict";

const gff_version = "r24";

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

/****************
 * FASTX reader *
 ****************/

class Fastx {
	constructor(fn) {
		this._file = new File(fn);
		this._last = 0;
		this._line = new Bytes();
		this._finished = false;
		this.s = new Bytes();
		this.q = new Bytes();
		this.n = new Bytes();
		this.c = new Bytes();
	}
	read() {
		var c, f = this._file, line = this._line;
		if (this._last == 0) { // then jump to the next header line
			while ((c = f.read()) != -1 && c != 62 && c != 64);
			if (c == -1) return -1; // end of file
			this._last = c;
		} // else: the first header char has been read in the previous call
		this.c.length = this.s.length = this.q.length = 0;
		if ((c = f.readline(this.n, 0)) < 0) return -1; // normal exit: EOF
		if (c != 10) f.readline(this.c); // read FASTA/Q comment
		if (this.s.capacity == 0) this.s.capacity = 256;
		while ((c = f.read()) != -1 && c != 62 && c != 43 && c != 64) {
			if (c == 10) continue; // skip empty lines
			this.s.set(c);
			f.readline(this.s, 2, this.s.length); // read the rest of the line
		}
		if (c == 62 || c == 64) this._last = c; // the first header char has been read
		if (c != 43) return this.s.length; // FASTA
		this.q.capacity = this.s.capacity;
		c = f.readline(this._line); // skip the rest of '+' line
		if (c < 0) return -2; // error: no quality string
		var size = this.s.length;
		while (f.readline(this.q, 2, this.q.length) >= 0 && this.q.length < size);
		f._last = 0; // we have not come to the next header line
		if (this.q.length != size) return -2; // error: qual string is of a different length
		return size;
	}
	close() {
		this._file.close();
		this.s.destroy();
		this.q.destroy();
		this.n.destroy();
		this.c.destroy();
	}
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

class Exon {
	constructor(est, een) {
		this.st = est, this.en = een;
	}
}

class Transcript {
	constructor(tid, ctg, strand) {
		this.tid = tid, this.ctg = ctg, this.strand = strand;
		this.st = this.en = this.cds_st = this.cds_en = -1;
		this.disp_name = this.target_name = null;
		this.type = null;
		this.gid = null, this.gname = null;
		this.pri = null;
		this.gff_flag = 0; // 1=MANE_Select, 2=Ensembl_canonical, 4=GENCODE_Primary, 8=readthrough_transcript
		this.score = -1;
		this.err = 0, this.done = false;
		this.exon = [];
		this.gff_cds = []; // for GFF/GTF parsing only; don't use elsewhere!
	}
	#update_gff_cds() { // compute cds_st and cds_en
		if (this.gff_cds.length == 0) return; // no CDS
		this.gff_cds = this.gff_cds.sort(function(a,b) {return a.st - b.st});
		const e0 = this.gff_cds[0];
		const e1 = this.gff_cds[this.gff_cds.length - 1];
		if (this.strand == "+") {
			this.cds_st = e0.phase < 0? e0.st : e0.st + e0.phase;
			this.cds_en = e1.phase < 0? e1.en : e1.st + e1.phase + Math.floor((e1.en - (e1.st + e1.phase)) / 3) * 3;
		} else if (this.strand == "-") {
			this.cds_st = e0.phase < 0? e0.st : e0.en - e0.phase - Math.floor(((e0.en - e0.phase) - e0.st) / 3) * 3;
			this.cds_en = e1.phase < 0? e1.en : e1.en - e1.phase;
		}
	}
	finish() { // FIXME: GTF and GFF3 are different on stop codon
		this.#update_gff_cds();
		if (this.exon.length == 0 && this.gff_cds.length != 0) // use CDS if exon is missing
			this.exon = this.gff_cds;
		this.gff_cds = []; // clear out
		if (this.exon.length == 0) return false;
		const a = this.exon.sort(function(a,b) {return a.st - b.st}); // sort by start
		this.exon = a;
		const st = a[0].st;
		const en = a[a.length - 1].en;
		if (this.st < 0) this.st = st;
		if (this.en < 0) this.en = en;
		if (st != this.st || en != this.en) this.err |= 1;
		for (let i = 1; i < a.length; ++i)
			if (a[i].st < a[i-1].en)
				this.err |= 2; // ERROR: overlapping exons
		if ((this.cds_st >= 0 && this.cds_st < st) || (this.cds_en >= 0 && this.cds_en > en))
			this.err |= 4; // ERROR: inconsistent cds_st/cds_en
		if (!this.disp_name)
			this.disp_name = this.gname && this.type? [this.tid, this.type, this.gname].join("|") : this.tid;
		if (this.err == 0) this.done = true;
		if (!this.done) this.perror();
		return this.done;
	}
	has_cds() {
		return this.cds_st >= 0 && this.cds_en > 0 && this.cds_en > this.cds_st? true : false;
	}
	cds_len() {
		if (!this.has_cds()) return 0;
		let len = 0;
		for (let i = 0; i < this.exon.length; ++i) {
			const e = this.exon[i];
			let max_st = e.st > this.cds_st? e.st : this.cds_st;
			let min_en = e.en < this.cds_en? e.en : this.cds_en;
			if (max_st < min_en) len += min_en - max_st;
		}
		return len;
	}
	trans_len() {
		let len = 0;
		for (let i = 0; i < this.exon.length; ++i)
			len += this.exon[i].en - this.exon[i].st;
		return len;
	}
	cut_to_cds() {
		if (!this.has_cds()) return false;
		let exon = [];
		for (let i = 0; i < this.exon.length; ++i) {
			const e = this.exon[i];
			let max_st = e.st > this.cds_st? e.st : this.cds_st;
			let min_en = e.en < this.cds_en? e.en : this.cds_en;
			if (max_st < min_en) exon.push(new Exon(max_st, min_en));
		}
		this.st = this.cds_st, this.en = this.cds_en;
		this.exon = exon;
		return true;
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

function* gff_read(fn) {
	const re_fmt = /^(##PAF\t\S+(\t\d+){3}\t[+-])|^(\S+(\t\d+){3}\t[+-])|^(\S+\t\d+\t\d+\t\S+\t\S+\t[+-])|^((\S+\t){3}\d+\t\d+\t\S+\t[+-])|^(\S+\t\d+\t\S+\t\d+\t\d+\t(\d+[MIDNSH=X])+)/;
	const re_gff = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|tag|Parent|Target)( "([^"]+)"|=([^;]+))/g;
	const re_cigar = /(\d+)([MIDNSHP=XFGUV])/g;

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
				v.exon.push(new Exon(v.st + off, v.st + off + len));
			}
			if (v.finish()) yield v;
			v = null;
		} else if (fmt == 2) { // GTF or GFF3
			let t = line.split("\t");
			if (t[2] != "exon" && t[2] != "CDS" && t[2] != "cds") continue;
			let m, tid = null, gid = null, type = "", gname = null, biotype = "", gff_flag = 0, pid = null, target = null;
			while ((m = re_gff.exec(t[8])) != null) {
				const key = m[1], val = m[3]? m[3] : m[4];
				if (key == "transcript_id") tid = val;
				else if (key == "transcript_type") type = val;
				else if (key == "transcript_biotype" || key == "gbkey") biotype = val;
				else if (key == "gene_name") gname = val;
				else if (key == "gene_id") gid = val;
				else if (key == "Parent") pid = val;
				else if (key == "Target") target = val.split(" ")[0];
				else if (key == "tag") {
					if (val == "MANE_Select") gff_flag |= 1;
					else if (val == "Ensembl_canonical") gff_flag |= 2;
					else if (val == "GENCODE_Primary") gff_flag |= 4;
					else if (val == "readthrough_transcript") gff_flag |= 8;
				}
			}
			if (tid == null) tid = pid;
			if (tid == null) continue;
			if (gname == null) gname = gid? gid : "*"; // infer gene name
			if (gid == null) gid = gname; // if gene_id missing, use gene name to identify a gene
			if (type == "" && biotype != "") type = biotype; // infer transcript type
			if (v == null || v.tid != tid) {
				if (v != null && v.finish()) yield v;
				v = new Transcript(tid, t[0], t[6]);
				v.type = type, v.gid = gid, v.gname = gname, v.target_name = target, v.gff_flag = gff_flag;
			}
			const st = parseInt(t[3]) - 1;
			const en = parseInt(t[4]);
			if (t[2] === "CDS" || t[2] === "cds") {
				const phase = /^\d+$/.test(t[7])? parseInt(t[7]) % 3 : -1;
				v.gff_cds.push({ st:st, en:en, phase:phase });
			} else if (t[2] === "exon") {
				v.exon.push(new Exon(st, en));
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
						v.exon.push(new Exon(st + e[i][0], st + e[i][1]));
				} else if (t[4] == '-') { // For protein-to-genome alignment, the coordinates are on the query strand. Need to flip them.
					const glen = en - st;
					for (let i = e.length - 1; i >= 0; --i)
						v.exon.push(new Exon(st + (glen - e[i][1]), st + (glen - e[i][0])));
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

function gff_select(g) // priority: MANE_select > Ensembl_canonical > GENCODE_Primary > longest_CDS > longest_transcript
{
	if (g.length == 1) return g[0];
	let a = [];
	for (let i = 0; i < g.length; ++i)
		a.push({ v:g[i], cds_len: g[i].cds_len(), trans_len: g[i].trans_len() });
	a = a.sort(function(x,y) { return x.cds_len != y.cds_len? x.cds_len - y.cds_len : x.trans_len - y.trans_len; });
	let sel_mane = null, sel_ens = null, sel_pri = null;
	for (let i = 0; i < a.length; ++i) {
		if (a[i].v.gff_flag & 1) sel_mane = a[i].v;
		if (a[i].v.gff_flag & 2) sel_ens = a[i].v;
		if (a[i].v.gff_flag & 4) sel_pri = a[i].v;
	}
	if (sel_mane) return sel_mane;
	else if (sel_ens) return sel_ens;
	else if (sel_pri) return sel_pri;
	else return a[0].v;
}

function* gff_read_select(fn)
{
	let g = [];
	for (let v of gff_read(fn)) {
		if (g.length > 0 && g[0].gid != v.gid) {
			yield gff_select(g);
			g = [];
		}
		g.push(v);
	}
	if (g.length > 0)
		yield gff_select(g);
}

function gff_cmd_all2bed(args)
{
	let pri_only = false, cds_only = false, disp_target_name = false, print_junc = false, print_ss = false, select1 = false, no_through = false;
	for (const o of getopt(args, "aptjs1r", [])) {
		if (o.opt == "-p") pri_only = true;
		else if (o.opt == "-a") cds_only = true;
		else if (o.opt == "-t") disp_target_name = true;
		else if (o.opt == "-j") print_junc = true;
		else if (o.opt == "-s") print_ss = true;
		else if (o.opt == "-1") select1 = true;
		else if (o.opt == "-r") no_through = true;
	}
	if (args.length == 0) {
		print("Usage: minigff.js all2bed [options] <in.file>");
		print("Options:");
		print("  -a       CDS only");
		print("  -1       one transcript per gene (GFF; only if clustered by gene)");
		print("  -p       only include primary alignments");
		print("  -j       print junctions/introns");
		print("  -s       print 3bp at splice sites");
		print("  -t       display Target name in BED");
		return;
	}
	let reader = select1? gff_read_select : gff_read;
	for (let v of reader(args[0])) {
		if (pri_only && !v.pri) continue;
		if (no_through && (v.gff_flag&8)) continue;
		if (cds_only) {
			if (!v.has_cds()) continue;
			v.cut_to_cds();
		}
		if (disp_target_name && v.target_name != null)
			v.disp_name = v.target_name;
		if (print_junc) {
			for (let i = 1; i < v.exon.length; ++i)
				print(v.ctg, v.exon[i-1].en, v.exon[i].st, ".", ".", v.strand);
		} else if (print_ss) {
			for (let i = 1; i < v.exon.length; ++i) {
				const st = v.exon[i-1].en, en = v.exon[i].st;
				if (v.strand == "+") {
					print(v.ctg, st, st+3, "D", ".", "+");
					print(v.ctg, en-3, en, "A", ".", "+");
				} else if (v.strand == "-") {
					print(v.ctg, en-3, en, "D", ".", "-");
					print(v.ctg, st, st+3, "A", ".", "-");
				}
			}
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
	add(v) { // add a Transcript object
		if (this.ctg[v.ctg] == null) this.ctg[v.ctg] = { exon:[], junc:[] };
		let p = this.ctg[v.ctg];
		for (let i = 0; i < v.exon.length; ++i)
			p.exon.push({ st:v.exon[i].st, en:v.exon[i].en });
		for (let i = 1; i < v.exon.length; ++i)
			p.junc.push({ st:v.exon[i-1].en, en:v.exon[i].st });
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
	const re_chr = /^(chr)?(\d{1,3}|[XYZW]|\d{1,3}[A-Z]|[IVX]+)$/;
	let print_all = false, print_err = false, first_only = false, chr_only = false, skip_last = false, skip_first = false, cds_only = false, eval_base = true;
	for (const o of getopt(args, "e1ctfdapB", [])) {
		if (o.opt == "-e") print_err = true;
		else if (o.opt == "-p") print_all = true;
		else if (o.opt == "-1") first_only = true;
		else if (o.opt == "-c") chr_only = true;
		else if (o.opt == "-f") skip_first = true;
		else if (o.opt == "-t") skip_first = skip_last = true;
		else if (o.opt == "-B") eval_base = false;
		else if (o.opt == "-d" || o.opt == "-a") cds_only = true;
	}
	if (args.length < 2) {
		print("Usage: minigff.js eval [options] <base.file> <test.file>");
		print("Options:");
		print("  -a      CDS only");
		print("  -1      only evaluate the first alignment of each TEST");
		print("  -c      only consider contigs matching /^(chr)?([0-9]{1,3}|[XYZW]|[0-9]{1,3}[A-Z])$/");
		print("  -f      skip the first exon in TEST for exon evaluation");
		print("  -t      skip the last exon in TEST for exon evaluation");
		print("  -B      skip base evaluation (less memory)");
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
	print("CC\tBN  nAnnoBase  nJointBase  sensitivity");
	print("CC\tBP  nAlnBase   nJointBase  specificity");
	print("CC");

	// load base annotation
	let base = new BaseIndex();
	for (let v of gff_read(args[0])) {
		if (chr_only && !re_chr.test(v.ctg)) continue;
		if (cds_only) {
			if (!v.has_cds()) continue;
			v.cut_to_cds();
		}
		base.add(v);
	}
	base.index();

	// evaluate test annotation
	let last_tid = null, tot_exon = 0, ann_exon = 0, nov_exon = 0, tot_junc = 0, ann_junc = 0, nov_junc = 0, n_multi = 0, n_test = 0;
	let qexon = {};
	for (let v of gff_read(args[1])) {
		if (chr_only && !re_chr.test(v.ctg)) continue;
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

/*********************
 * Extract sequences *
 *********************/

function gff_print_fasta(name, seq, line_len)
{
	print(`>${name}`);
	if (line_len <= 0) {
		print(seq);
	} else {
		const s = seq.toString();
		for (let i = 0; i < s.length; i += line_len)
			print(s.substr(i, line_len));
	}
}

function gff_cmd_getseq(args)
{
	let cds_only = false, select1 = false, to_translate = false, line_len = 80, skip_flawed_aa = false;
	for (const o of getopt(args, "a1tl:f", [])) {
		if (o.opt == "-a") cds_only = true;
		else if (o.opt == "-1") select1 = true;
		else if (o.opt == "-t") cds_only = to_translate = true;
		else if (o.opt == "-l") line_len = parseInt(o.arg);
		else if (o.opt == "-f") skip_flawed_aa = true;
	}
	if (args.length < 2) {
		print("Usage: minigff.js getseq [options] <anno.bed> <seq.fa>");
		print("Options:");
		print("  -a      only extract CDS");
		print("  -1      one transcript per gene (GFF; only if clustered by gene)");
		print("  -t      translate CDS (forcing -a)");
		print("  -f      skip protein sequences with stop codons");
		print(`  -l INT  line length [${line_len}]`);
		return 1;
	}
	let reader = select1? gff_read_select : gff_read;
	let bed = {};
	for (let v of reader(args[0])) {
		if (cds_only) {
			if (!v.has_cds()) continue;
			v.cut_to_cds();
		}
		if (bed[v.ctg] == null) bed[v.ctg] = [];
		bed[v.ctg].push(v);
	}

	// prepare translation table
	const codon2aa = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX";
	const stop_code = '*'.charCodeAt(0);
	let codon2aa_table = new Uint8Array(new ArrayBuffer(65));
	for (let i = 0; i <= 64; ++i)
		codon2aa_table[i] = codon2aa.charCodeAt(i);

	let nt4_table = new Uint8Array(new ArrayBuffer(256));
	for (let i = 0; i < 256; ++i) nt4_table[i] = 4;
	nt4_table['A'.charCodeAt(0)] = nt4_table['a'.charCodeAt(0)] = 0;
	nt4_table['C'.charCodeAt(0)] = nt4_table['c'.charCodeAt(0)] = 1;
	nt4_table['G'.charCodeAt(0)] = nt4_table['g'.charCodeAt(0)] = 2;
	nt4_table['T'.charCodeAt(0)] = nt4_table['t'.charCodeAt(0)] = 3;

	let len, fx = new Fastx(args[1]);
	let buf = new Bytes();
	while ((len = fx.read()) >= 0) {
		const seq = new Uint8Array(fx.s.buffer);
		const ctg = fx.n.toString();
		if (bed[ctg] == null) continue;
		const v = bed[ctg];
		for (let i = 0; i < v.length; ++i) {
			buf.length = 0;
			for (let j = 0; j < v[i].exon.length; ++j)
				buf.set(seq.slice(v[i].exon[j].st, v[i].exon[j].en).buffer);
			if (v[i].strand == "-") k8_revcomp(buf);
			const disp_name = v[i].disp_name? v[i].disp_name : v[i].tid;
			if (to_translate) { // translate to protein sequence
				const s = new Uint8Array(buf.buffer);
				let n_stop = 0;
				if (s.length % 3 != 0)
					warn(`Warning: ${disp_name} - CDS length is not a multiple of 3; continue anyway`);
				let aa_seq = new Bytes(Math.floor((buf.length + 2) / 3));
				let aa_arr = new Uint8Array(aa_seq.buffer);
				for (let j = 0; j < aa_seq.length; ++j) {
					const c1 = nt4_table[s[j*3]], c2 = nt4_table[s[j*3+1]], c3 = nt4_table[s[j*3+2]];
					const codon = c1 < 4 && c2 < 4 && c3 < 4? (c1<<4|c2<<2|c3) : 64;
					const aa = codon2aa_table[codon];
					if (aa == stop_code) {
						if (j == aa_seq.length - 1) // if the last codon is a stop codon, ignore it
							break;
						++n_stop;
					}
					aa_arr[j] = aa;
				}
				if (n_stop > 0) warn(`Warning: ${disp_name} contains ${n_stop} stop codon(s)`);
				if (!skip_flawed_aa || n_stop == 0)
					gff_print_fasta(disp_name, aa_seq, line_len);
				aa_seq.destroy();
			} else { // print nucleotide sequence
				gff_print_fasta(disp_name, buf, line_len);
			}
		}
	}
	buf.destroy();
	fx.close();
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
		print("  getseq         extract transcript sequences");
		print("  version        print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == "all2bed") gff_cmd_all2bed(args);
	else if (cmd == "eval") gff_cmd_eval(args);
	else if (cmd == "getseq") gff_cmd_getseq(args);
	else if (cmd == "version") {
		print(gff_version);
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
