//
// PREPARE INTERVALS
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
include { BUILD_INTERVALS                                        } from '../../../modules/local/build_intervals/main'
include { CREATE_INTERVALS_BED                                   } from '../../../modules/local/create_intervals_bed/main'
include { GATK4_INTERVALLISTTOBED                                } from '../../../modules/nf-core/gatk4/intervallisttobed/main'

workflow PREPARE_INTERVALS {
    take:
    fasta_fai    // mandatory [ fasta_fai ]
    intervals    // [ params.intervals ]

    main:
    versions = Channel.empty()

    if (!intervals) {
        BUILD_INTERVALS(fasta_fai.map{it -> [ [ id:it.baseName ], it ] })

        intervals_combined = BUILD_INTERVALS.out.bed

        CREATE_INTERVALS_BED(intervals_combined.map{ meta, path -> path }).bed

        intervals_bed = CREATE_INTERVALS_BED.out.bed

        versions = versions.mix(BUILD_INTERVALS.out.versions)
        versions = versions.mix(CREATE_INTERVALS_BED.out.versions)
    } else {
        intervals_combined = Channel.fromPath(file(intervals)).map{it -> [ [ id:it.baseName ], it ] }
        intervals_bed = CREATE_INTERVALS_BED(file(intervals)).bed

        versions = versions.mix(CREATE_INTERVALS_BED.out.versions)

        // If interval file is not provided as .bed, but e.g. as .interval_list then convert to BED format
        if (intervals.endsWith(".interval_list")) {
            GATK4_INTERVALLISTTOBED(intervals_combined)
            intervals_combined = GATK4_INTERVALLISTTOBED.out.bed
            versions = versions.mix(GATK4_INTERVALLISTTOBED.out.versions)
        }
    }

    // Now for the intervals.bed the following operations are done:
    // 1. Intervals file is split up into multiple bed files for scatter/gather
    // 2. Each bed file is indexed

    // 1. Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals
    intervals_bed = intervals_bed.flatten()
        .map{ intervalFile ->
            def duration = 0.0
            for (line in intervalFile.readLines()) {
                final fields = line.split('\t')
                if (fields.size() >= 5) duration += fields[4].toFloat()
                else {
                    start = fields[1].toInteger()
                    end = fields[2].toInteger()
                    duration += (end - start) / params.nucleotides_per_second
                }
            }
            [ duration, intervalFile ]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
        .flatten().collate(2).map{ duration, intervalFile -> intervalFile }.collect()
        // Adding number of intervals as elements
        .map{ it -> [ it, it.size() ] }
        .transpose()

    intervals_bed_combined        = intervals_combined.map{meta, bed -> bed }.collect()

    emit:
    // Intervals split for parallel execution
    intervals_bed                 // [ intervals.bed, num_intervals ]
    // All intervals in one file
    intervals_bed_combined        // [ intervals.bed ]
    versions               // [ versions.yml ]
}
