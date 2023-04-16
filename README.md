# RNAseq_BALLGOWN_PIPELINE

To use this script the following command with command line argmuments indicating file names must be incldued:
```ruby
sh pipeline.sh <referemce_sample.fna> <sample_base_name> <name_of_sample_gtf.gtf>
```
The usage of this script requieres modifications before being ran on the shell terminal. These additional requirements are: 
1. Files: This script must be edited where any <> are found. Within the script anytime a <> is found, it countains a description of what it should represent, once modifieed remove the <>. 

2. Directory: This script must be ran within a directory containing : 
- All fastq.gz files
- Reference species/sample transcript for alignment = first command line argument 
- General feature format or general transfer format for mapping. GTF format is the one requireed but this script has a command to modify from GFF > GTF. This script also includes a command to consistently format GFF files downloaded from internet reference. Both are commented out inittially.  


3. Package dependency. This script requieres the following modules to be loaded: 
- fastqc
- hisat2
- samtools
- stringtie
- gffread (optional)



### STEP 1: Run fastqc
FASTQC is done as a quality check on the RNA reads. If the quality is low, it requieres an additional step of trimming - not displayed in this pipeline - 

```ruby
module load fastqc
mkdir ./FASTQC
fastqc *.fastq.gz  -o ./FASTQC  
```



For cleaning up FNA files downloaded from internet, and ensuring they only contain identifier in header row which is neccesary for creating indices, uncomment the command: 

```ruby
awk -F "." '{print $1}' "$fna_file" > GCFedited.fna
```
genomic.fna is the reference_sample.fna file the user provides as a command line argument


### STEP 2: Run Hisat2--build 
This step is an alignment of your reads with genome indexes. This pipeline first creates the index based on a DNA fasta file and then aligns reads to the created index
Usage:
- hisat2-build [options]* <reference_in> <ht2_base>
- <reference_in> A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if -c is specified, the sequences themselves
- <ht2-base>  The basename of the index files to write. By default, hisat2-build writes files named NAME.1.ht2, NAME.2.ht2 where NAME is <ht2_base>
-- large-index  Force hisat2-build to build a large index, even if the reference is less than ~ 4 billion nucleotides long.
- f The reference input files (specified as <reference_in>) are FASTA files (usually having extension .fa, .mfa, .fna or similar).


```ruby 
module load hisat2/2.0.5
hisat2-build --large-index  -f GCFedited.fna <ht2-base>
       
```

### Step 3: Run hisat2 alignment step 
This step contains <> that should be changed to the last string pattern in common between all samples before file denotation (i.e. .fna .gff .gtf etc ) 

HISAT2 main usage alignment
hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> } [-S <hit>]
- dta
       is for downstream aplications such as Stringtie
-p is for processors being used
-S is for the output sam file
-x path for the indices built and being used for alignment
-1&2 these are the paired end reads to be aligned 
```ruby
 for fq1 in ./*<R1_001>.fastq.gz  
    do

    base=$(basename $fq1 <R1_001>.fastq.gz)


    fq1=./${base}<R1_001>.fastq.gz
    fq2=./${base}<R2_001>.fastq.gz

hisat2 -p 8 --dta -x ./ -1 $fq1 -2 $fq2 -S ${base}.${sample_base_name}.sam
done

```


### Step 4: converting SAM files to BAM files
To do anything meaningful with alignment data you must swith from SAM to its binary counterpart BAM file. This binary format is much easier for computer programs such as StringTie to work with.
Basic usage: 
$ samtools <command> [options] 
- sort this sorts our bam sam files into bam files for more applications 
-o gives the output file name

       
```ruby
module load samtools
for i in ./*.sam  
    do
    base=$(basename $i .sam)  
        

    samtools sort  ${i}  -o ${base}.bam
done
```

In case the user cannot fine a gft file fore the following step, this command can convert the gff file into a gft file 
```ruby
module load gffread/
gffread genomic.gff -T -o genomic.gtf
```


### STEP 5 Run featureCounts - StringTie
StringTie assembles RNAseq alignments from HiSAT into potential transcripts. It can use a gtf file as a guide for the assembly. 
An assembled gtf file will be generated for each sample by using StringTie. 
The command :
 
```ruby 
ls *_ST_.gtf > mergelist.txt  
```
Is a wildcard to take all gtf files made from StringTie and include them all within the mergelist.txt file
The last command within this step takes the recently created merged .txt file and returns a proper merged gtf file called stringtie_merged.gtf

StringTie Main usage:
stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
- p number of proccesors
- G reference genome 
- o output gtf file
- l label, this is the prefix for the name of the output transcripts. Defualt: STRG 

      
module load stringtie/1.3.3
```ruby
module load stringtie
ffor file in *.bam
do
       tag=${file%.bam}
stringtie -p 8 -G "$genomic_gtf" -o $tag_ST_.gtf -l ${sample_base_name}  $tag.bam
done

ls *ST_.gtf > mergelist.txt 

           
stringtie --merge -p 8 -G "$genomic.gtf" -o stringtie_merged.gtf mergelist.txt

```


### STEP 6 Generating count table for ballgown:
To use ballgown, a count table from all the gtf files must be made. Stringtie will compare each sample agianst merged assembly to espimate transcript abundance
```ruby

module load stringtie

for file in *.bam;

    do tag=${file%.bam};
    mkdir ./ballgown/$tag/$tag.gtf
    stringtie -e -B -p 8 -G stringtie_merged.gtf -o ./ballgown/$tag/$tag.gtf $tag.bam

done

```
Now you will have neccesary output files to use in ballgown analysis, such files are : 
1. e_data.ctab
2. e2t.ctab
3. i_data.ctab
4. i2t.ctab
5. t_data.ctab
You should have these 5 files within a directory called ballgown and within subdirectories based of the sample names 
# RNAseq_ballgown
