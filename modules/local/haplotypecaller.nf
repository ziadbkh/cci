
process HAPLOTYPECALLER {
    tag "${sample_id} - ${interval.getName()}"
    
    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    
    input:
    tuple val(sample_id), path(bam), path(bam_bai)
    path reference_genome
    each path (interval) 
    
    output:
    tuple val(sample_id), path("${sample_id}*.g.vcf.gz"), emit: gvcf
    tuple val(sample_id), path("${sample_id}*.vcf.gz.tbi"), emit: gvcf_tbi

    script:
    gvcf_basename = "${sample_id}.${interval.getName()}.g.vcf.gz"
    refernce_genome_main = reference_genome.find { it.getName().endsWith(".fa") }

    """
    gatk --java-options "-Xmx8000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=2" \
        HaplotypeCaller \
        -R ${refernce_genome_main} \
        -I ${bam} \
        -L ${interval} \
        -O ${gvcf_basename} \
        -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
        -stand-call-conf 15.0 \
        -GQB 5 -GQB 10 -GQB 15 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 \
        -ERC GVCF

    """
}
