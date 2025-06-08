## Getting Started
```sh
# download
git clone https://github.com/lh3/minigff
cd minigff

# obtain the k8 binary, or install via bioconda
wget https://zenodo.org/records/11357203/files/k8-1.2.tar.bz2
tar -jxvf k8-1.2.tar.bz2
cp k8-1.2/k8-`uname -m`-`uname -s` $HOME/bin/k8  # assuming $HOME/bin in $PATH

# conversion
./minigff.js all2bed test/gc47.gtf.gz > out.bed     # convert to BED12
./minigff.js all2bed -1 test/gc47.gtf.gz > out.bed  # one transcript per gene

# evaluation
./minigff.js eval -a1 test/gc47.gtf.gz test/mp.paf.gz
```

## Introduction

minigff is a script that parses gene annotation formats (GTF, GFF3 and BED12)
and spliced alignment formats (SAM and PAF), extracts information and compares
transcript structures. It seamlessly reads multiple formats, though for
GTF/GFF3, transcripts are required to be grouped by genes. minigff is writen in
a dialect of Javascript and depends on the [k8][k8] Javascript engine. Precompiled
k8 binaries can downloaded [from zenodo][k8-dl] or installed [via bioconda][k8-bc].

minigff grew out of my needs for processing and evaluating spliced alignment.
It unifies and *deprecates* `gff2bed.js` in [minisplice][msp] and the
`splice2bed`, `gff2bed`, `gff2junc`, `junceval` and `exoneval` subcommands in
[paftools.js][paftools]. [gffread][gffr] and [gffcompare][gffc] are similar to
minigff in functionality and provide more features. They are more widely used.

[msp]: https://github.com/lh3/minisplice
[paftools]: https://github.com/lh3/minimap2/tree/master/misc
[k8]: https://github.com/attractivechaos/k8
[k8-dl]: https://zenodo.org/records/11357203
[gffr]: https://github.com/gpertea/gffread
[gffc]: https://github.com/gpertea/gffcompare
[k8-bc]: https://bioconda.github.io/recipes/k8/README.html
