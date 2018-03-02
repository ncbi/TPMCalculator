TPMCalculator
===

This program calculates the TPM values for the genes, transcripts, exons and introns 
from NGS RNA-Seq aligned reads (BAM files). It also 
compile a set of files per samples with the raw counts, TPM for exons and 
introns individually for the gene and transcripts.    

This program parse a GTF file and create a list of genes overlapping the exon 
regions among the transcripts and creating pure intronic regions. After this, 
the genes are compared with the other in order to overlap exons and introns. 
The final construction include exonic regions and pure intronic regions for 
each gene.

After that, the program parses RNA-Seq aligned reads (BAM files) in a directory 
or a single file allocating the reads into the genomic features for the before 
created genes and transcripts.
 
Once all reads for a BAM file are allocated the program compute the TPM value 
for the exons, introns and for the union of the exons (gene and transcript) and 
the union of the introns (intronic region for genes and transcripts). 

## Output files

### Output file: _gene.out
This file is created per sample and includes the TPM values at gene or 
transcript level. 

The columns are:

Exon_Length\tExon_Count_Reads\tExon_TPM\tIntron_Length\tIntron_Count_Reads\tIntron_TPM"

1. Gene_Id
2. Chr
3. Length (Transcript length)
4. Count_Reads (Reads assigned to the gene or transcript)
5. TPM (TPM value for the gene or transcript)
6. Length of non overlapping gene exons
7. Reads assigned to the non overlapping gene exons
8. TPM for the non overlapping gene exons
9. Length of non overlapping gene introns
10. Reads assigned to the non overlapping gene introns
11. TPM for the non overlapping gene introns
12. Exon_Length (Sum of all exonic regions)
13. Exon_Count_Reads (All reads assigned to the exonic regions)
14. Exon_TPM (TPM for all the exonic regions)
15. Intron_Length (Sum of all intronic regions)
16. Intron_Count_Reads (All reads assigned to the intronic regions)
17. Intron_TPM (TPM for all the intronic regions)

### Output file: _transcript.out
This file is created per sample and includes the TPM values at gene or 
transcript level. 

The columns are:

1. Gene_Id
2. Transcript_Id (Not present in the gene file)
3. Chr
4. Length (Transcript length)
5. Count_Reads (Reads assigned to the gene or transcript)
6. TPM (TPM value for the gene or transcript)
7. Exon_Length (Sum of all exonic regions)
8. Exon_Count_Reads (All reads assigned to the exonic regions)
9. Exon_TPM (TPM for all the exonic regions)
10. Intron_Length (Sum of all intronic regions)
11. Intron_Count_Reads (All reads assigned to the intronic regions)
12. Intron_TPM (TPM for all the intronic regions)

### Output file: _[gene|transcript].ent
This file is per sample and includes the TPM values for each exon and intron
at gene or transcript level.
 
The columns are:

1. Gene_Id
2. Transcript_Id (Not present in the gene file)
3. Chr
4. Type ( feature type: exon or intron)
5. Type_Number (consecutive number starting from 1)
6. start (feature start coordinate)
7. end (feature end coordinate)
8. Length (feature length)
9. Count_Reads (reads assigned to the feature)
10. TPM (TPM calculated for the feature)

### Output file: [genes|transcripts]_data_per_samples.txt
This file is created when multiple BAM files are processed from a directory.
It includes the TPM values for exonic region and intronic region 
calculated for each gene or transcript for each sample. 

The columns are (using the same names than before):

1. Gene_Id
2. Transcript_Id (Not present in the gene file)
3. sample1_exon
4. sample1_intron
5. sample2_exon
6. sample2_intron

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

Usage: ./bin/TPMCalculator -g GTF_file [-d BAM_files_directory|-i BAM_file] 

    ./bin/TPMCalculator options:

    -v    Print info
    -h    Display this usage information.
    -g    GTF file
    -d    Directory with the BAM files
    -b    BAM file
    -k    Gene key to use from GTF file. Default: gene_id
    -t    Transcript key to use from GTF file. Default: transcript_id
    -c    Smaller size allowed for an intron created for genes. Default: 16. We recommend to use the reads length
    -p    Use only properly paired reads. Default: No. Recommended for paired-end reads.

## Credits

Roberto Vera Alvarez, PhD

Emails: veraalva@ncbi.nlm.nih.gov, r78v10a07@gmail.com

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

