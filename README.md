# scafSeqContigInfo

### Contents
#### scafSeqContigInfo.py
Python script that calculates assembly statistics at the scaffold level as well as a summary set of statistics.  
#### scafSeqContigInfo.sh
wrapper shell script for scafSeqContigInfo.py with specific options to calculate assembly summary statistics, skipping statistics on each individual scaffold. It calculates the statistics on all scaffolds and after excluding scaffolds smaller than three length cutoffs (201 nt, 300 nt, 1,000 nt).  

### Usage
$ git clone https://github.com/calacademy-research/scafSeqContigInfo.git  
$ cd scafSeqContigInfo  
\# Run "scafSeqContigInfo.sh" on a file in FASTA format  
$ ./scafSeqContigInfo.sh \<assembly.fasta\> \<results_file_prefix\>  

### scafSeqContigInfo.py options
$ python scafSeqContigInfo.py  
```
Usage: scafSeqContigInfo.py <filename.scafSeq> [ [-x] [-s <smallest_scaff_to_include>] | -All ] [-1 .. -99]

       Looks through scaffolds in the scafSeq file (by default those 1K or larger) and displays overview info,
       By default reports one line for scaffold, summarizing scaffold len, contig lens, num contigs, largest contig, N50, %Ns
       the output is sorted by the scaffold length largest to smallest (this length is the contig bases + N gaps).

       -s option: by default only scaffold lengths >= 1K are shown, this can be changed to as low as 201 by -s option
       -S option: same as -s and list of Longest Contigs will include Scaffold names
       -x option: this excludes the individual info lines for each scaffold and only shows the overview info lines
       -nc option: do not show longest and shortest contig info lines
       -All option: shows overview info for all scaffolds and singleton contigs, individual info lines are not shown
       -1 .. -99, e.g. -4, sets gap threshold to number of N's allowed in contig before splitting it (0 by default)
```
### scafSeqContigInfo.sh description
$ ./scafSeqContigInfo.sh  
```
scafSeqContigInfo.sh does runs of the script to compute scaffold info using several scaffold size cutoffs.
It takes the scaffold fasta file and an assembly prefix for the output file as arguments, e.g.

      scafSeqContigInfo.sh *scafSeq asm28
         # this example would output to a file named asm28_scafContigInfo.txt
```

### Citing

#### Authorship

Code author: James B. Henderson, jhenderson@calacademy.org  
README.md authors: Zachary R. Hanna, James B. Henderson  

#### Version 1.0.0
[![DOI](https://zenodo.org/badge/67550118.svg)](https://zenodo.org/badge/latestdoi/67550118)
