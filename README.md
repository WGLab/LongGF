# LongGF
A computational algorithm and software tool for fast and accurate detection of gene fusion by long-read transcriptome sequencing

## Background
Long-read RNA-Seq techniques can generate reads that encompass a large proportion or the entire mRNA or cDNA molecules, so they are expected to address inherited limitations of short-read RNA-Seq techniques that generate only 50-150bp reads. However, there is a general lack of software tools for gene fusion detection from long-read RNA-seq data, which takes into account of the higher error rates and the possible alignment errors for long-read data. Here, we proposed a fast computational tool, LongGF, to efficiently detect candidate gene fusion from long-read RNA-seq data, including cDNA sequencing data and direct mRNA sequencing data. 

## Input
A bam file of long-read transcriptome sequencing data sorted by name and a GTF file for the definition of genes. **Please note that you bam must be sorted by NAME but NOT by POSITION.**

## Requirements
You need to have C++ (GCC > 4.8.5) and HTSlib to compile and run the program. 

## Installation
There are several steps before you can run the program.
1. `git clone https://github.com/WGLab/LongGF`
2. `cd LongGF/bin`
3. `g++ -g -O3 -std=c++11 -I ./include -L ./lib -Wl,--enable-new-dtags,-rpath,"\$ORIGIN"/lib -lhts -o LongGF _com_fun_.c _gtf_struct_.c get_gfFrombam.c -Wl,--no-as-needed`

Then, you will have `LongGF/bin/LongGF` to run.

## Usage
You can run `./LongGF` to get help documents about the parameters belw:
```
Usage: ./LongGF <input_bam> <input_gtf> <min-overlap-len> <bin_size> <min-map-len> [pseudogene:0(default)/1] [Secondary_alignment:0(default)] [min_sup_read:2(default)]
```
where
```
<input_bam>: A bam file. Compulsory.
<input_gtf>: A GTF file. Compulsory. 
<min-overlap-len>: The minimum length of an alignment record overlap with a gene. Compulsory. Recommend: 100.
<bin_size>: The bin size during discretization. Compulsory. Recommended: 30, 50 or 100.
<min-map-len>: A minimum length of an alignment record against the reference genome. Compulsory. Recommend: 100, 200 or 300.
[pseudogene:0(default)/1]: Optional. Default=0: Not use pseudogene from the GTF file.
[Secondary_alignment:0(default)]: Optional. Default=0: not use secondary alignment.
[min_sup_read:2(default)]: Optional. Default=2. 
```

## Contact
If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/WGLab/LongGF/issues). They would also be helpful to other users.

## Reference
Qian Liu, Yu Hu, Andres Stucky, Li Fang, Jiang F. Zhong, Kai Wang. LongGF: computational algorithm and software tool for fast and accurate detection of gene fusion by long-read transcriptome sequencing. Submitted.


