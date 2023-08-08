/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCreatepanelofnormals.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.dict,
    params.fasta,
    params.fasta_fai,
    params.intervals,
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

if(params.tools && params.tools.split(',').contains("cnvkit")){
    if(!params.sequencing_method){
        error("--sequencing_method is required for CNVKit. See the CNVKit documentation for details.")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CNVKIT_BATCH                      } from '../modules/nf-core/cnvkit/batch/main'
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../modules/nf-core/gatk4/createsomaticpanelofnormals/main'
include { GATK4_GENOMICSDBIMPORT            } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_MERGEVCFS                   } from '../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MUTECT2                     } from '../modules/nf-core/gatk4/mutect2/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { PREPARE_INTERVALS                 } from '../subworkflows/local/prepare_intervals/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CREATEPANELOFNORMALS {

    ch_versions = Channel.empty()

    input         = Channel.fromSamplesheet("input")
    fasta         = params.fasta     ? Channel.fromPath(params.fasta).first()     : Channel.empty()
    fai           = params.fasta_fai ? Channel.fromPath(params.fasta_fai).first() : Channel.empty()
    dict          = params.dict      ? Channel.fromPath(params.dict).first()      : Channel.empty()
    intervals_all = params.intervals ? Channel.fromPath(params.intervals).first() : Channel.empty()

    fasta = fasta.map{ it -> [[id:it.baseName], it]}
    fai   = fai.map{ it -> [[id:it.baseName], it]}
    dict  = dict.map{ it -> [[id:it.baseName], it]}

    if(params.tools && params.tools.split(',').contains('mutect2')){

        PREPARE_INTERVALS(fai, params.intervals)

        intervals = PREPARE_INTERVALS.out.intervals_bed

        // Combine input and intervals for spread and gather strategy
        input_intervals = input.combine(intervals)
            // Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
            .map{ meta, input_list, input_index_list, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input_list, input_index_list, intervals ] }

        GATK4_MUTECT2(input_intervals,
                    fasta,
                    fai,
                    dict,
                    [],[],[],[])

        ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

        // Figuring out if there is one or more vcf(s) from the same sample
        vcf_branch = GATK4_MUTECT2.out.vcf.branch{
            // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

        // Figuring out if there is one or more tbi(s) from the same sample
        tbi_branch = GATK4_MUTECT2.out.tbi.branch{
            // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

        vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }
                                        .groupTuple()
        GATK4_MERGEVCFS(vcf_to_merge, dict)
        ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

        vcf = Channel.empty().mix(GATK4_MERGEVCFS.out.vcf, vcf_branch.no_intervals)
        tbi = Channel.empty().mix(GATK4_MERGEVCFS.out.tbi, tbi_branch.no_intervals)

        vcf_joint = vcf.map{ meta, vcf -> [[id: 'joint'], vcf]}.groupTuple()
        tbi_joint = tbi.map{ meta, tbi -> [[id: 'joint'], tbi]}.groupTuple()

        ch_genomicsdb_input = vcf_joint.join(tbi_joint).combine(intervals_all)
                                        .map{ meta, vcf, tbi, intervals ->
                                            [meta, vcf, tbi, intervals, [], []]
                                        }

        GATK4_GENOMICSDBIMPORT(ch_genomicsdb_input,
                                [],
                                [],
                                [])
        ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

        GATK4_CREATESOMATICPANELOFNORMALS(GATK4_GENOMICSDBIMPORT.out.genomicsdb,
                                            fasta,
                                            fai,
                                            dict)
        ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions)
    }

    if(params.tools && params.tools.split(',').contains('cnvkit')){

        pooled_normal = input.map{meta, cram, crai -> [[id:"normal"], cram]}
                            .groupTuple()
        pooled_normal.view()

        CNVKIT_BATCH(pooled_normal, fasta, intervals_all)

        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

    }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCreatepanelofnormals.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCreatepanelofnormals.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
