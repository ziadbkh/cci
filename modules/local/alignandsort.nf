process ALIGNANDSORT {
    tag "${sample_id}"    
    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'zhangb1/biobambam2-samtools-picard-bwa-samblaster-sambamba' }"

    
    input:
    tuple val(sample_id), val(id), val(paired_end), path(fastqs)
    path (reference_genome)
    
    output:
    tuple val(sample_id), path("${output_bam_basename}.bam"), emit: aligned_and_sorted_bam
    tuple val(sample_id), path("${output_bam_basename}.bam.bai"), emit: aligned_and_sorted_bam_bai
    
    script:
    output_bam_basename = "${sample_id}.${id}.dedup.sorted"
    metrics_filename="${sample_id}.${id}.bam.dupmetrics"
    
    fastqs_str = fastqs.join(" ")
    bamsorThreads = Math.min(12,task.cpus)
    refernce_genome_main = reference_genome.find { it.getName().endsWith(".fa") }

    """
    bwa mem \
        -t ${task.cpus} \
        -R @RG\\\\tID:"${sample_id}"\\\\tLB:None\\\\tPL:ILLUMINA\\\\tPU:NA\\\\tSM:"${sample_id}" \
        -v 3 -Y \
        ${refernce_genome_main} \
        ${fastqs_str} | \
    bamsormadup \
        threads=${bamsorThreads} \
        inputformat=sam \
        outputformat=bam \
        optminpixeldif=2500 \
        reference=${refernce_genome_main} \
        tmpfile=tmp_${output_bam_basename} \
        fragmergepar=${bamsorThreads} \
        M=${metrics_filename} \
        indexfilename=${output_bam_basename}.bam.bai \
        > ${output_bam_basename}.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        bamsormadup: \$(echo \$(bamsormadup --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
    END_VERSIONS
    """

}