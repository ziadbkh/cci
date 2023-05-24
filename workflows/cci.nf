/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCci.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { ALIGNANDSORT                  } from '../modules/local/alignandsort'
include { MERGEBAMS                     } from '../modules/local/mergebams.nf'
include { HAPLOTYPECALLER               } from '../modules/local/haplotypecaller.nf'
include { MERGEVCFS                     } from '../modules/local/mergevcfs.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow CCI {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    
    //INPUT_CHECK (    ch_input )
    //ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    


    /*Channel
        .fromFilePairs(params.input + '*_R{1,2}*.fastq.gz', checkIfExists:true)
        .map{[it[0], it[1][0].simpleName, it[1]]}
        .set{ch_input_fastq}
    */
    
    Channel
        .fromPath( params.fasta, checkIfExists:true).first().set{ch_reference_genome}

    Channel
        .fromPath( params.fasta + '.{bwt,sa,ann,amb,pac}', checkIfExists:true).toSortedList().set{ch_reference_genome_extra_bwa}
    
    
    ch_reference_genome.map{"${it.getParent()}/${it.getBaseName()}.dict"}
    .combine(Channel
            .fromPath( params.fasta + '.fai', checkIfExists:true)
    )
    .set{ch_reference_genome_extra_gatk}


    Channel
        .fromPath( params.intervals + '*scattered.interval_list' , checkIfExists:true).set{ch_intervals}

    Channel
        .fromPath( params.input, checkIfExists:true )
        .splitCsv(header: true)
        .map { 
            row -> {
            if (row.fastq2)
                {    
                    [row.sample_id, row.id, row.type.toLowerCase(), true, [row.fastq1, row.fastq2] ] 
                }else{
                    [row.sample_id, row.id, row.type.toLowerCase(), false, [row.fastq1] ]
                }
            }
           
        }
        .set{ch_meta_complete}
    
    ch_meta_complete
    .map{[it[0], it[1], it[3], it[4]]}
    .set{ch_meta}
    
    
    ch_meta_complete
    .filter{it[2] != "t" && it[2] != "tumor" && it[2] != "0"}
    .map{[it[0], it[2]]}
    .set{ch_normal_samples}

    ALIGNANDSORT(
        ch_meta,
        ch_reference_genome.combine(ch_reference_genome_extra_bwa).flatten().toSortedList()
    )
    
    ALIGNANDSORT.out.aligned_and_sorted_bam
    .groupTuple().join(
        ALIGNANDSORT.out.aligned_and_sorted_bam_bai.groupTuple()
    )
    .branch {
        singelton: it[1].size() <= 1
        multiple: it[1].size() > 1
    }
    .set {
        ch_sample_bams
    }
    
    
    MERGEBAMS(
        ch_sample_bams.multiple
    )
    
    
    MERGEBAMS.out.merged_bam.join(MERGEBAMS.out.merged_bam_bai)
    .mix(
        ch_sample_bams
        .singelton
        .map{[it[0], it[1][0], it[2][0]]}
    ).set{
        ch_all_bams
    }
    
    ch_all_bams
    .join(ch_normal_samples)
    .map{[it[0], it[1]]}
    .set{ch_bams}

    HAPLOTYPECALLER(
        ch_bams,
        ch_reference_genome
            .combine(ch_reference_genome_extra_gatk)
            .flatten()
            .toSortedList(),
        ch_intervals

    )
    
    HAPLOTYPECALLER.out.gvcf
    .groupTuple()
    .join(
        HAPLOTYPECALLER.out.gvcf_tbi
        .groupTuple()
    ).set{ch_sample_gvcf}
    
    MERGEVCFS(
        ch_sample_gvcf
    )
    
    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    */
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
