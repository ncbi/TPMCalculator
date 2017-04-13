TPMCalculator
===

This program calculates the TPM values for the exons and introns. It also 
compile a set of files per samples with the raw counts, TPM for exons and 
introns individually and for the whole transcript.    

The program parses all BAM files in a directory and uses a GTF file to allocate
the reads into the genomic features. Once all reads for a BAM file are allocated
the program compute the TPM value for the exons, introns and for the union of 
the exons (transcript) and the union of the introns (intronic region). 

Further analysis is executed to compute the mean of the TPM values between two
conditions. The program includes two option -c (control group) and 
-t (treated group) that are the prefix pattern used to identify the control 
and treated groups from the BAM file names. 

For example, let's suppose that the conditions prefix patterns are NBC for 
control samples and CLL for treated samples. 

The BAM files can be named as:
    NBC_1072821_accepted_hits.sorted.bam
    NBC_1072822_accepted_hits.sorted.bam
    CLL_1072723_accepted_hits.sorted.bam
    CLL_1072724_accepted_hits.sorted.bam

The program will calculate a CLL_exon TPM and CLL_intron TPM from the mean of 
all CLL exon and intron TPM values respectively. It will do the same with the 
NBC samples.   

## Output files

### Output file: .out
This file is created per sample and includes the TPM values at transcript level. 

The columns are:

1. Gene_Id
2. Transcript_Id
3. Chr
4. Length (Transcript length)
5. Count_Reads (Reads assigned to the transcript)
6. TPM (TPM value for transcript)
7. Exon_Length (Sum of all exonic regions)
8. Exon_Count_Reads (All reads assigned to the exonic regions)
9. Exon_TPM (TPM for all the exonic regions)
10. Intron_Length (Sum of all intronic regions)
11. Intron_Count_Reads (All reads assigned to the intronic regions)
12. Intron_TPM (TPM for all the intronic regions)

### Output file: .ent
This file is per sample and includes the TPM values for each exon and intron.
 
The columns are:

1. Gene_Id
2. Transcript_Id
3. Chr
4. Type ( feature type: exon or intron)
5. Type_Number (consecutive number starting from 1)
6. start (feature start coordinate)
7. end (feature end coordinate)
8. Length (feature length)
9. Count_Reads (reads assigned to the feature)
10. TPM (TPM calculated for the feature)

### Output file: transcript_all_per_sample.txt
This file include the TPM values for exonic region and intronic region 
calculated for each transcript for each sample. 

The columns are (using the same names than before):

1. Gene_Id
2. Transcript_Id
3. NBC_1072821_accepted_hits.sorted_exon
4. NBC_1072821_accepted_hits.sorted_intron
5. NBC_1072822_accepted_hits.sorted_exon
6. NBC_1072822_accepted_hits.sorted_intron
7. CLL_1072723_accepted_hits.sorted_exon
8. CLL_1072723_accepted_hits.sorted_intron
9. CLL_1072724_accepted_hits.sorted_exon
10. CLL_1072724_accepted_hits.sorted_intron

### Output file: transcript_count_per_samples.txt
This file include the raw count values for the whole transcript (exonic regions
only) calculated for each transcript for each sample. 

The columns are (using the same names than before):

1. Gene_Id
2. Transcript_Id
5. NBC_1072821_accepted_hits.sorted
6. NBC_1072822_accepted_hits.sorted
7. CLL_1072723_accepted_hits.sorted
8. CLL_1072724_accepted_hits.sorted

### Output file: intron_count_per_samples.txt
This file include the raw count values for the introns calculated for each 
transcript for each sample. 

The columns are (using the same names than before):

1. Gene_Id
2. Transcript_Id
3. Intron_Number (consecutive number starting from 1)
4. Intron_Id (transcript id concatenated with intron number. 
   Transcript id: NR_028057, Intorn 1 Id: NR_028057_1)
5. NBC_1072821_accepted_hits.sorted
6. NBC_1072822_accepted_hits.sorted
7. CLL_1072723_accepted_hits.sorted
8. CLL_1072724_accepted_hits.sorted

## Requirements

### BAMTools

Clone the BAMTools repository from GitHub: https://github.com/pezmaster31/bamtools

Compile it on this way and set the environment variables for TPMCalculator:

    cd bamtools
    mkdir build
    cd build
    cmake ..
    make
    cd ..
    export BAMTOOLS_DIR=`pwd`
    export CPPFLAGS="-I $BAMTOOLS_DIR/include"
    export LDFLAGS="-L $BAMTOOLS_DIR/lib -Wl,-rpath,$BAMTOOLS_DIR/lib"

That's it. BAMTools was compiled and the env variables were set for compiling
TPMCalculator.

## Installation

After the instalation of BAMTools go to the TPMCalculator folder and do make:

    make

A bin folder will be created with the TPMCalculator executable.

## Usage

Usage: ./bin/TPMCalculator -g GTF_file -d BAM_files_directory -c NBC -t CLL

    ./bin/TPMCalculator options:

    -v    Print info
    -h    Display this usage information.
    -g    GTF file
    -d    Directory with the BAM files
    -c    Control identification string. This is a string pattern in the control BAM file names
    -t    Treated identification string. This is a string pattern in the treated BAM file names

## Credits

Roberto Vera Alvarez, PhD

Emails: veraalva@ncbi.nlm.nih.gov, r78v10a07@gmail.com

===
# Public Domain notice

## National Center for Biotechnology Information.

This software is a "United States Government Work" under the terms of the United States
Copyright Act. It was written as part of the authors' official duties as United States
Government employees and thus cannot be copyrighted. This software is freely available
to the public for use. The National Library of Medicine and the U.S. Government have not
placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy and reliability
of the software and data, the NLM and the U.S. Government do not and cannot warrant the
performance or results that may be obtained by using this software or data. The NLM and
the U.S. Government disclaim all warranties, express or implied, including warranties
of performance, merchantability or fitness for any particular purpose.

Please cite NCBI in any work or product based on this material.

