![FixAME](/images/FixAME.png)


**FixAME** is a program for finding and fixing local assembly errors, as well as identifying whether these errors are possible chimeras.

## Installation

First of all, you need to make sure all dependecies are installed at your $PATH

[samtools](https://github.com/samtools/samtools) v.1.12+ <br>
[bcftools](https://github.com/samtools/bcftools) v.1.12+<br> 
[BBMAP](https://sourceforge.net/projects/bbmap/) -- specially `bbmap.sh` and `filterbyname.sh`

-- Python 3.6+ and the dependecies

[pysam](https://github.com/pysam-developers/pysam) v.0.15.4+<br>
[xopen](https://pypi.org/project/xopen/) v.0.9.0 <br> 
[pandas](https://pypi.org/project/pandas/) v.1.0.3+ <br>
[Biopython](https://biopython.org/wiki/Download) v.1.77

After that you can clone this repository to the desired location

```
git clone https://github.com/LiviaMoura/FixAME.git
```

## Running FixAME

All you need to run FixAME is a sample or a folder with bins with extension `.fa`|`.fasta`|`.fna` and a pair-ended reads (R1 and R2)

FixAME ignores fasta sequences smaller than `1000bp` by default (800bp minimum). Theses sequences won`t be extended or have the errors located.

FixAME **CAN'T** deal with `bins` named using the same name differing with .1, .2, .3|.a, .b, .c| etc (ex.: bin.1, bin.2, bin.3)

### Processing a single genome or a metagenome fasta

```
FixAME.py --fasta sample.fasta -r1 R1.fastq -r2 R2.fastq --output_dir THIS_FOLDER 
````

### Processing a folder containing bins

```
FixAME.py --bins THIS_FOLDER_CONTAINS_BINS -r1 R1.fastq -r2 R2.fastq --output_dir THIS_FOLDER 
```

### Error finder
If you don't want to run the fixing on the found errors, you can only run the error_finder step. This is return only the local assembly errors that FixAME could identify.

```
FixAME.py error_finder -f sample.fa -r1 R1.fastq -r2 R2.fastq -o TEST_FOLDER
```










