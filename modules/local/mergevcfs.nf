
process MERGEVCFS {
    tag "${sample_id}"
    
    conda (params.enable_conda ?  "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
        
    
    input:
    tuple val(sample_id), path(gvcfs), path(gvcf_tbis)

    output:
    tuple val(sample_id), path("${output_vcf_name}"), emit: merged_vcf
    tuple val(sample_id), path("${output_vcf_name}.tbi"), emit: merged_vcf_tbi
    
    script:
    output_vcf_name = "${sample_id}.hc.merged.hs38.g.vcf.gz"
    
    """
    gatk --java-options "-Xmx8000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=2" \
        MergeVcfs \
        ${gvcfs.collect{ "--INPUT $it"}.join(' ')} \
        --OUTPUT ${output_vcf_name}
    """
}