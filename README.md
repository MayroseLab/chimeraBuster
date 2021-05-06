# chimeraBuster
chimeraBuster is a tool for correcting chimeric gene annotations.  
Please **do not** use yet - still under development.  

## Installation and running
chimeraBuster was tested on Linux, but should also run on Windows/Mac.  
### Requirements
* Git
* python (any version)
* conda (Anaconda/Miniconda, compatible with the python version)
### Installation
```
$ git clone git@github.com:MayroseLab/chimeraBuster.git
$ cd chimeraBuster/
$ conda env create -f env.yml
$ conda activate chimeraBuster
$ python chimeraBuster.py -h
```
### Running a test data set
```
$ python chimeraBuster.py -g test/genes.gff -f test/genome.fasta -t test/transcripts.fasta -o test_result
```
