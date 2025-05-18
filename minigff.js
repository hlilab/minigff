#!/usr/bin/env k8

"use strict";

const gff_version = "r2";

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

function iit_sort_copy(a) {
	a.sort((x, y) => (x.st - y.st));
	const b = [];
	for (let i = 0; i < a.length; ++i)
		b.push({ st: a[i].st, en: a[i].en, max: 0, data: a[i].data });
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

/**************
 *** gfform ***
 **************/

class Transcript {
	constructor(tid, ctg, strand) {
		this.tid = tid, this.ctg = ctg, this.strand = strand;
		this.st = this.en = this.cds_st = this.cds_en = -1;
		this.disp_name = null;
		this.type = null;
		this.gid = null, this.gname = null;
		this.err = 0, this.done = false;
		this.exon = [];
	}
	finish() {
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
		if (this.cds_st < 0) this.cds_st = st;
		if (this.cds_en < 0) this.cds_en = en;
		if (this.cds_st < st || this.cds_en > en) this.err |= 4;
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
		return [this.ctg, this.st, this.en, this.disp_name, ".", this.strand, this.cds_st, this.cds_en, ".", this.exon.length, lens.join(",")+",", offs.join(",")+","];
	}
}

function* gff_read(fn) {
	const re_fmt = /^(##PAF\t\S+(\t\d+){3}[+-])|^(\S+(\t\d+){3}[+-])|^(\S+\t\d+\t\d+\t\S+\t\S+\t[+-])|^((\S+\t){3}\d+\t\d+\t\S+\t[+-])/;
	const re_gff = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|tag)( "([^"]+)"|=([^;]+));/g;
	let v = null;
	for (const line of k8_readline(fn)) {
		let m;
		if ((m = re_fmt.exec(line)) == null) continue;
		let ft = m[5]? 1 : m[6]? 2 : m[3]? 3 : m[1]? 4 : -1; // 1=BED, 2=GFF, 3=PAF, 4=##PAF
		if (ft < 0) continue;
		if (ft == 1) {
			let t = line.split("\t");
			v = new Transcript(t[3], t[0], t[5]);
			v.st = parseInt(t[1]), v.en = parseInt(t[2]);
			v.cds_st = parseInt(t[6]), v.cds_en = parseInt(t[7]);
			const n_exon = parseInt(t[9]);
			const lens = t[10].split(",", n_exon);
			const offs = t[11].split(",", n_exon);
			for (let i = 0; i < n_exon; ++i) {
				const off = parseInt(offs[i]), len = parseInt(lens[i]);
				v.exon.push({ st: v.st + off, en: v.st + off + len });
			}
			if (v.finish()) yield v;
			v = null;
		} else if (ft == 2) {
			let t = line.split("\t");
			if (t[2] != "exon" && t[2] != "CDS") continue;
			let m, tid = null, gid = null, type = "", gname = null, biotype = "", ens_canonical = false;
			while ((m = re_gff.exec(t[8])) != null) {
				const key = m[1], val = m[3]? m[3] : m[4];
				if (key == "transcript_id") tid = val;
				else if (key == "transcript_type") type = val;
				else if (key == "transcript_biotype" || key == "gbkey") biotype = val;
				else if (key == "gene_name") gname = val;
				else if (key == "gene_id") gid = val;
				else if (key == "tag" && val == "Ensembl_canonical") ens_canonical = true;
			}
			if (tid == null) continue;
			if (gname == null) gname = gid? gid : "*"; // infer gene name
			if (gid == null) gid = gname; // if gene_id missing, use gene name to identify a gene
			if (type == "" && biotype != "") type = biotype; // infer transcript type
			if (v == null || v.tid != tid) {
				if (v != null && v.finish()) yield v;
				v = new Transcript(tid, t[0], t[6]);
				v.type = type, v.gid = gid, v.gname = gname;
			}
			const st = parseInt(t[3]) - 1;
			const en = parseInt(t[4]);
			if (t[2] === "CDS") {
				v.cds_st = v.cds_st >= 0 && v.cds_st < st? v.cds_st : st;
				v.cds_en = v.cds_en >= 0 && v.cds_en > en? v.cds_en : en;
			} else if (t[2] === "exon") {
				v.exon.push({ st:st, en:en });
			}
		}
	}
	if (v != null && v.finish()) yield v;
}

function gff_cmd_read(args)
{
	for (const o of getopt(args, "", [])) {
	}
	if (args.length == 0) {
		print("Usage: minigff.js read [options] <in.file>");
		return;
	}
	for (let v of gff_read(args[0])) {
		print(v.bed().join("\t"));
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
		print("  read           parse BED12/GFF/GTF/PAF");
		print("  version        print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'read') gff_cmd_read(args);
	else if (cmd == 'version') {
		print(gff_version);
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
