/*
 * bismark subworkflow
 */
include { ARIOC_ALIGN                               } from '../../modules/nf-core/AriocP/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED } from '../../modules/nf-core/samtools/sort/main'
include { ARIOC_DEDUPLICATE                         } from '../../modules/nf-core/bismark/deduplicate/main'
include { ARIOC_METHYLATIONEXTRACTOR                } from '../../modules/nf-core/bismark/methylationextractor/main'
include { ARIOC_COVERAGE2CYTOSINE                   } from '../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                              } from '../../modules/nf-core/bismark/report/main'
include { ARIOC_SUMMARY                             } from '../../modules/nf-core/bismark/summary/main'

workflow ARIOC {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    arioc_index      // channel: /path/to/BismarkIndex/
    skip_deduplication // boolean: whether to deduplicate alignments
    cytosine_report    // boolean: whether the run coverage2cytosine

    main:
    versions = Channel.empty()


    /*
     * Align with arioc
     */
    ARIOC_ALIGN (
        reads,
        arioc_index
    )
    versions = versions.mix(ARIOC_ALIGN.out.versions)

    /*
     * Sort raw output BAM
     */
    SAMTOOLS_SORT_ALIGNED(
        ARIOC_ALIGN.out.bam,
    )
    versions = versions.mix(SAMTOOLS_SORT_ALIGNED.out.versions)

    if (skip_deduplication) {
        alignments = ARIOC_ALIGN.out.bam
        alignment_reports = ARIOC_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    } else {
        /*
        * Run deduplicate_bismark
        */
        ARIOC_DEDUPLICATE( ARIOC_ALIGN.out.bam )

        alignments = ARIOC_DEDUPLICATE.out.bam
        alignment_reports = ARIOC_ALIGN.out.report.join(ARIOC_DEDUPLICATE.out.report)
        versions = versions.mix(ARIOC_DEDUPLICATE.out.versions)
    }

    /*
     * Run bismark_methylation_extractor
     */
    ARIOC_METHYLATIONEXTRACTOR (
        alignments,
        arioc_index
    )
    versions = versions.mix(ARIOC_METHYLATIONEXTRACTOR.out.versions)


    /*
     * Run coverage2cytosine
     */
    if (cytosine_report) {
        ARIOC_COVERAGE2CYTOSINE (
            ARIOC_METHYLATIONEXTRACTOR.out.coverage,
            arioc_index
        )
        versions = versions.mix(ARIOC_COVERAGE2CYTOSINE.out.versions)
    }

    /*
     * Generate bismark sample reports
     */
    ARIOC_REPORT (
        alignment_reports
            .join(ARIOC_METHYLATIONEXTRACTOR.out.report)
            .join(ARIOC_METHYLATIONEXTRACTOR.out.mbias)
    )
    versions = versions.mix(ARIOC_REPORT.out.versions)

    /*
     * Generate bismark summary report
     */
    ARIOC_SUMMARY (
        ARIOC_ALIGN.out.bam.collect{ it[1].name }.ifEmpty([]),
        alignment_reports.collect{ it[1] }.ifEmpty([]),
        alignment_reports.collect{ it[2] }.ifEmpty([]),
        ARIOC_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
        ARIOC_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
    )
    versions = versions.mix(ARIOC_SUMMARY.out.versions)

    /*
     * MODULE: Run samtools sort
     */
    SAMTOOLS_SORT_DEDUPLICATED (
        alignments
    )
    versions = versions.mix(SAMTOOLS_SORT_DEDUPLICATED.out.versions)

    /*
     * Collect MultiQC inputs
     */
    ARIOC_SUMMARY.out.summary.ifEmpty([])
        .mix(alignment_reports.collect{ it[1] })
        .mix(alignment_reports.collect{ it[2] })
        .mix(ARIOC_METHYLATIONEXTRACTOR.out.report.collect{ it[1] })
        .mix(ARIOC_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] })
        .mix(ARIOC_REPORT.out.report.collect{ it[1] })
        .set{ multiqc_files }

    emit:
    bam        = SAMTOOLS_SORT_ALIGNED.out.bam        // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup      = SAMTOOLS_SORT_DEDUPLICATED.out.bam   // channel: [ val(meta), [ bam ] ] ## sorted, possibly deduplicated BAM
    mqc        = multiqc_files                        // path: *{html,txt}
    versions                                       // path: *.version.txt
}
