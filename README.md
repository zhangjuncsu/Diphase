# Diphase: phase diploid genome using heterozygous variations
Diphase is a phasing tool designed for diploid genome using heterozyous variations. Diphase can only work on the primary/alternate assembly format.
## dependencies
Diphase is compiled with C++11. We compiled Diphase using [gcc-12.2](https://gcc.gnu.org/gcc-12/) on CentOS 7. Python3 is needed to run the pipeline.
- [minimap2](https://github.com/lh3/minimap2/)
- [samtools](https://github.com/samtools/samtools)
- [BWA-MEM](https://github.com/lh3/bwa)
- [Clair3](https://github.com/HKU-BAL/Clair3)
### dependencies for htslib
- zlib
- libbz2
- liblzma
- libcurl
- libcrypto
## Installation instructions
### Install Diphase from bioconda (recommended)
```
conda create -n <env_name> diphase -c bioconda
```
### Install Diphase from GitHub
Download the latest code from GitHub:
```
git clone https://github.com/zhangjuncsu/Diphase.git
cd Diphase/src & make
cd ..
export PATH=`pwd`/bin:$PATH
```
Diphase can be found in ./bin and the python scripts can be found in ./script.
## Usage
```
python /pathto/pipeline.py --pri <primary assembly> --alt <alternate assembly> --rdfname <reads> --hic1 <Hi-C mate-pair 1> --hic2 <Hi-C mate-pair 2> --model <clair3 model path> --type [clr | hifi | ont] -d <out directory> -t <threads>
```
### Options
| | | |
| :--- | :--- | :---|
| --pri | \<FILE> | primary assembly file name 
--alt | \<FILE> | alternate assembly file name |
--rdfname | \<FILE> | raw reads file name |
--hic1 | \<FILE> | Hi-C pair 1 file name |
--hic2 | \<FILE> | Hi-C pair 2 file name |
--model | \<FILE> | model path used to call SNP by Clair3 |
--type | \<STR> | sequencing technology [ont] |
--seed | \<INT> | seed for random function default 10000 |
--iter | \<INT> | iteration for phasing algorithm default 1000 |
-d | \<DIR> | directory to save the output files |
-t | \<INT> | number of threads default 1 |
-q | \<INT> | mapping quality to filter Hi-C mapping default 1 |
-p | \<STR> | prefix of the output files default "phasing" |

Diphase can only work on the primary/alternate asembly format, now. We are implementing the code to make Diphase work on the other assembly format, such as dual assembly format.
