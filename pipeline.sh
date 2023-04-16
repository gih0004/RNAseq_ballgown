#!/bin/bash





#To use this script the following command with command line argmuments indicating file names must be incldued:
# $sh pipeline.sh <referemce_sample.fna> <sample_base_name> <name_of_sample_gtf.gtf>

#The usage of this script requieres modifications before being ran on the shell terminal. These additional requirements are: 
#1. Files: This script must be edited where any <> are found. Within the script anytime a <> is found, it countains a description of what it should represent, once modifieed remove the <>. 

#2. Directory: This script must be ran within a directory containing : 
#a. all fastq.gz files
#b. reference species/sample transcript for alignment = first command line argument 
#c. General feature format or general transfer format for mapping. GTF format is the one requireed but this script has a command to modify from GFF > GTF. This script also includes
#a command to consistently format GFF files downloaded from internet reference.  


#3. Package dependency. This script requieres the following modules to be loaded: 
#a. fastqc
#b. hisat2
#c. samtools
#d. stringtie
#e. gffread (optional)

fna_file="$1"
sample_base_name="$2"
genomic_gtf="$3" 

module load fastqc/0.11.9
module load hisat2/2.2.1
module load stringtie/2.1.6
module load samtools/1.15 
#module load gffread/0.9.8 


if [ "$#" -ne 3 ]; then
    echo "Error: Expected three command line arguments." > that_was_not_right.txt  # Creates error file if command line arguments are not three
    echo "Aborting the program."

#Step 1 RUN FASTQC 

#The flag -o indicated the output directory 
mkdir ./FASTQC
fastqc *.fastq.gz  -o ./FASTQC  



#STEP 2: RUN HISAT2 hisat2-build 

#<_R1_001> and <_R2_001> must be changed to the last common string between samples before file format denotation

#Usage:
#hisat2-build [options]* <reference_in> <ht2_base>
#<reference_in> A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if -c is specified, the sequences themselves
#<ht2-base> The basename of the index files to write. By default, hisat2-build writes files named NAME.1.ht2, NAME.2.ht2 where NAME is <ht2_base>
#--large-index  Force hisat2-build to build a large index, even if the reference is less than ~ 4 billion nucleotides long.
#-f The reference input files (specified as <reference_in>) are FASTA files (usually having extension .fa, .mfa, .fna or similar).


#awk -F "." '{print $1}' "$fna_file" > GCFedited.fna
mkdir ./${sample_base_name}_index 
cd ./${sample_base_name}_index             # this creates a directory specifically for the index
hisat2-build --large-index  -f ${fna_file} ${sample_base_name} 
 
cd ..             # this changes back to the main directory 

#STEP 3: Run hisat2 alignment step 

#HISAT2 main usage alignment
#hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> 
#-dta is for downstream aplications such as Stringtie
#-p is for processors being used
#-S is for the output sam file
#-x path for the indices built and being used for alignment


for fq1 in ./*<R1_001>.fastq.gz #<>This should be changed into whichever string you have last in sample names common between all samples 
    do

    base=$(basename $fq1 <R1_001>.fastq.gz)


    fq1=./${base}<R1_001>.fastq.gz
    fq2=./${base}<R2_001>.fastq.gz

    hisat2 -p 8 --dta -x ./${sample_base_name}_index -1 $fq1 -2 $fq2 -S ${base}.${sample_base_name}.sam
done



# STEP 4: RUN SAMTOOLS


#Basic usage: 
#$ samtools <command> [options] 

for i in ./*.sam  
    do
    base=$(basename $i .sam)  
        

    samtools sort  ${i}  -o ${base}.bam
done


#STEP 5: Run featureCounts- String Tie 
#When using StringTie, a GTF file instead of a GFF file is requiered. These can be quickly converted following the following commands: 
#$module load gffread/0.9.8                                 # Or equivelant within HPC in use 
#$gffread <genomic/file_name>.gff -T -o genomic.gtf


#The command :
# $ls *_ST_.gtf > mergelist.txt  # is a wildcard to take all gtf files made from StringTie and include them all within the mergelist.txt file

#The last command within this step takes the recently created merged .txt file and returns a proper merged gtf file called stringtie_merged.gtf

#StringTie Main usage:
#stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
#-p number of proccesors
#-G reference genome 
#-o output gtf file
#-l label, this is the prefix for the name of the output transcripts. Defualt: STRG 

for file in *.bam
    do
    tag=${file%.bam}
    stringtie -p 8 -G "$genomic_gtf" -o $tag_ST_.gtf -l ${sample_base_name}  $tag.bam
done

ls *ST_.gtf > mergelist.txt 

           
stringtie --merge -p 8 -G "$genomic.gtf" -o stringtie_merged.gtf mergelist.txt

"""
Step 5: Gene Count table :

To use ballgown, a count table from all the gtf files must be made. Stringtie will compare each sample agianst merged assembly to estimate transcript abundance
"""


for file in *.bam;

    do tag=${file%.bam};
    mkdir ./ballgown/$tag/$tag.gtf
    stringtie -e -B -p 8 -G stringtie_merged.gtf -o ./ballgown/$tag/$tag.gtf $tag.bam

done


