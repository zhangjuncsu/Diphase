# Diphase: phasing Nanopore genome assembly by integrating heterozygous variations and Hi-C data
Diphase is a phasing tool designed for diploid genome using heterozyous variations. Diphase can only work on the primary/alternate assembly format.
## dependencies
Diphase is compiled with C++11. We compiled Diphase using [gcc-12.2](https://gcc.gnu.org/gcc-12/) on CentOS 7. Python3 is needed to run the pipeline. We recommend users install the dependencies using conda.
- [minimap2](https://github.com/lh3/minimap2/)
- [samtools](https://github.com/samtools/samtools)
- [BWA-MEM](https://github.com/lh3/bwa)
- [Clair3](https://github.com/HKU-BAL/Clair3)

Note: script ```run_clair3.sh``` should be included in PATH.
### dependencies for htslib
- zlib
- libbz2
- liblzma
- libcurl
- libcrypto
## Installation instructions
### Install Diphase from bioconda (recommended)
```
conda create -n <env_name> -c bioconda -c conda-forge diphase clair3 minimap2=2.22
```
model can be found in /path/to/conda/envs/<env_name>/bin/models/
### Docker pre-built image
A pre-built docker image is available at [DockerHub](https://hub.docker.com/repository/docker/jun67/diphase/general)
```
docker run -it \
-v ${INPUT_DIR}:/mnt \
-v ${OUTPUT_DIR}:/mnt/${OUTPUT_DIR} \
jun67/diphase:latest \
pipeline.py \
phase \
--pri ${INPUT_DIR}/primary.fasta \      ## primary assembly file
--alt ${INPUT_DIR}/alternate.fasta \    ## alternate assembly file
--rdfname ${INPUT_DIR}/reads.fastq.gz \ ## raw reads file
--hic1 ${INPUT_DIR}/hic_R1.fastq.gz \   ## Hi-C reads 1 file
--hic2 ${INPUT_DIR}/hic_R2.fastq.gz \   ## Hi-C reads 2 file
--model /opt/models/ont \               ## model path used to call SNP by Clair3
--type ont \                            ## sequencing technology [ont]
-d ${OUTPUT_DIR} \                      ## directory to save the output files
-t ${THREADS} \                         ## maximal threads to be used
--dump_filtered                         ## save filtered Hi-C mapping
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
## Testing
Download the testing data from [Google Drive](https://drive.google.com/file/d/1rvvWr4t4ZjbuJPP6PrLujh6FHxRmKE5e/view?usp=drive_link). The example data are used for input. Then run the demo to test whether diphase has been successfully installed.
```
tar -zxf data.tar.gz
/path/to/Diphase/script/pipeline.py phase --pri /path/to/data/primary.fasta --alt /path/to/data/alternate.fasta --rdfname /path/to/data/subread.fastq.gz --hic1 /path/to/data/HiC1.fastq.gz --hic2 /path/to/data/HiC2.fastq.gz --model <clair3 model path> --type ont -d <out directory> -t <threads> --dump_filtered
```
[Here](https://drive.google.com/file/d/1KiybiVVkIzygzCfZL9bL69yGEiyISVKG/view?usp=drive_link) is the output for the example data. The phased results can be found in ```diphase/phase/phasing.result.txt```. The final phased assemblies ```phasing.hap1.fasta``` and ```phasing.hap2.fasta``` can be found in the directory ```diphase```.
## Usage
### prepare your datasets
Prepare the primary contigs file, alternate contigs file, raw reads file, and Hi-C reads file. All the files can be gzip compressed. Then run
```
/path/to/Diphase/script/pipeline.py phase --pri <primary assembly> --alt <alternate assembly> --rdfname <reads> --hic1 <Hi-C mate-pair 1> --hic2 <Hi-C mate-pair 2> --model <clair3 model path> --type [clr | hifi | ont] -d <out directory> -t <threads> --dump_filtered
```

- SNPs called by Clair3 are in the files ```<out directory>/clair3/clair3pri/merge_output.vcf.gz``` and ```<out directory>/clair3/clair3alt/merge_output.vcf.gz```.
- Filtered Hi-C alignments are in the file ```<out directory>/phase/phasing.hic.filtered.bam```.
- Detected switches are in the file ```<out directory>/phasing.fixed.switch```.
- The phasing results are in the file ```<out directory>/phase/phasing.result.txt```.
- The final phased assemblied are in the files ```<out directory>/${prefix}.hap1.fasta``` and ```<out directory>/${prefix}.hap2.fasta```.
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
--dump_filtered | | dump filtered Hi-C mapping

Diphase can only work on the primary/alternate asembly format, now. We are implementing the code to make Diphase work on the other assembly format, such as dual assembly format.
