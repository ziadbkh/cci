
process MERGEBAMS {
    tag "${sample_id}"

    conda (params.enable_conda ?  "bioconda::sambamba=1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0--h98b6b92_0':
        'quay.io/biocontainers/sambamba:1.0--h98b6b92_0' }"
        
    input:
    tuple val(sample_id), path (bam_list), path(bam_bai_list)
    
    output:
    tuple val(sample_id), path ("${output_bam}"), emit: merged_bam
    tuple val(sample_id), path ("${output_bam}.bai"), emit: merged_bam_bai
    
    script:
    output_bam="${sample_id}.dedup.merged.sorted.hs38.bam"
    
    """
    sambamba merge -t ${task.cpus} ${output_bam} ${bam_list.join(" ")}
    """
}