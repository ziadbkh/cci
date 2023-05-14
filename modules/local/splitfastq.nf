
process SplitFastq {
    
    conda "bioconda::fastqsplitter=1.2.0 "
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqsplitter%3A1.2.0--py39hbf8eff0_3' :
        'quay.io/biocontainers/fastqsplitter:1.2.0--py310h1425a21_3' }"
        

    input:
    each line from lines

    output: 
    tuple val(sampleID),fastq1,fastq2 into CreateInputFileChannel

    script:
    list = line.split('\t')
    sampleID = list[0]
    fastq1 = list[1]
    fastq2 = list[2]
    """
    singularity exec \${SIFDIR}/fastqsplitter.sif fastqsplitter -i ${fastq1} \
        -o ${sampleID}_1.R1.fastq.gz \
        -o ${sampleID}_2.R1.fastq.gz \
        -o ${sampleID}_3.R1.fastq.gz &
    singularity exec \${SIFDIR}/fastqsplitter.sif fastqsplitter -i ${fastq2} \
        -o ${sampleID}_1.R2.fastq.gz \
        -o ${sampleID}_2.R2.fastq.gz \
        -o ${sampleID}_3.R2.fastq.gz &
    wait
    """
}