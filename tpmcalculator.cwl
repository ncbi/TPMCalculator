#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
    DockerRequirement:
        dockerImageId: biocontainers/tpmcalculator:0.0.1
        dockerFile:
          $include: https://raw.githubusercontent.com/ncbi/TPMCalculator/master/Dockerfile

inputs:
    out_stdout:
        type: string
    out_stderr:
        type: string
    g:
        type: File
        inputBinding:
            position: 1
            prefix: -g
        doc: |
            GTF file
    d:
        type: Directory?
        inputBinding:
            position: 2
            prefix: -d
        doc: |
           Directory with the BAM files 
    b:
        type: File?
        inputBinding:
            position: 2
            prefix: -b
        doc: |
            BAM file
    k:
        type: string?
        inputBinding:
            position: 3
            prefix: -k
        doc: |
            Gene key to use from GTF file. Default: gene_id
    t:
        type: string?
        inputBinding:
            position: 3
            prefix: -t
        doc: |
            Transcript key to use from GTF file. Default: transcript_id
    c:
        type: int?
        inputBinding:
            position: 3
            prefix: -c
        doc: |
            Smaller size allowed for an intron created for genes. Default: 16. We recommend to use the reads length
    p:
        type: boolean?
        inputBinding:
            position: 3
            prefix: -p
        doc: |
            Use only properly paired reads. Default: No. Recommended for paired-end reads.
    q:
        type: int?
        inputBinding:
            position: 3
            prefix: -q
        doc: |
            Minimum MAPQ value to filter out reads. Default: 0. This value depends on the aligner MAPQ value.
    o:
        type: int?
        inputBinding:
            position: 3
            prefix: -o
        doc: |
            Minimum overlap between a reads and a feature. Default: 8.
    e:
        type: boolean?
        inputBinding:
            position: 3
            prefix: -e
        doc: |
           Extended output. This will include transcript level TPM values. Default: No. 

outputs:
    out_stdout:
        type: stdout
    out_stderr:
        type: stderr
    out_output:
        type: File[]
        outputBinding:
            glob: "*.out"
    ent_output:
        type: File[]
        outputBinding:
            glob: "*.ent"
    uni_output:
        type: File[]
        outputBinding:
            glob: "*.uni"

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["TPMCalculator"]

