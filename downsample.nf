/* ======================================================================== *
 *   Nextflow pipeline for downsampling BAMs                                *
 *                                                                          *
 * ======================================================================== */

help_message = "TODO"


/* ------------------------------------------------------------------------ *
 *   Config
 * ------------------------------------------------------------------------ */

// show help if --help flag is used
params.help = false
if ( params.help ) { log.info help_message; exit 0 }

// initialise variables
params.outdir      = false
params.input_bam   = false
params.num_repeats = false
params.num_reads   = false
params.range_reads = false

// check for required input data
if ( !params.input_bam )   { log.error "ERROR: input_bam is required"; exit 1 }
if ( !params.num_repeats ) { log.error "ERROR: num_repeats is required"; exit 1 }

// make sure exactly one of num_reads and range_reads have been inputted
if ( !params.range_reads && !params.num_reads ) {
    log.error "ERROR: one of num_reads or range_reads is required"; exit 1

} else if ( params.range_reads && params.num_reads ) {
    log.error "ERROR: only one of num_reads or range_reads is required, not both"; exit 1

// parse input
} else if ( params.range_reads && !params.num_reads ) { 
    range_min  = params.range_reads.split(',')[0].toInteger()
    range_max  = params.range_reads.split(',')[1].toInteger()
    range_step = params.range_reads.split(',')[2].toInteger()
    target_reads = (range_min..range_max).step(range_step)
    
} else if ( !params.range_reads && params.num_reads ) {
    target_reads = params.num_reads.split(',')
}

// TODO - check optional input data
//if ( !params.outdir )   { 
//   log.warn "output directory not set, defaulting to current working directory";
//    params.outdir = baseDir
//}

// load file objects
input_bam = file("$params.input_bam", checkIfExists: true)

// create summary for log file
def summary = [:]
summary["Input BAM"]   = params.input_bam
summary["Reads"]       = target_reads.join(',')
summary["Num repeats"] = params.num_repeats
summary["Output dir"]  = params.outdir

log.info """
==============================================================================
Inputs:
${ summary.collect { k, v -> "  ${k.padRight(12)}: $v"}.join("\n") }
==============================================================================
"""

/* ------------------------------------------------------------------------ *
 *   Pipeline
 * ------------------------------------------------------------------------ */

// define initial channels
// channel for each target depth
Channel
    .from( target_reads )
    .map{ "${it}reads" }
    .set{ depths_ch }


// channel for each repeat
Channel
    .from( 1..params.num_repeats.toInteger() )
    .map{ "rep$it" }
    .set{ reps_ch }


/* ------------------------------------------------------------------------ *
 *   Check input BAM
 * ------------------------------------------------------------------------ */

// TODO - limit to bed region


// count total reads in input
process get_input_read_count {
    conda "$baseDir/env/downsample.yml"

    input:
    file(bam) from input_bam

    output:
    env(READ_COUNT) into input_bam_read_count_ch

    script:
    """
    READ_COUNT=\$(samtools view -c $bam)
    """
}


// are there enough reads to meet num_reads?
// TODO? - check inputs are integers
process input_bam_qc {
    conda "$baseDir/env/downsample.yml"

    input:
    val(input_count) from input_bam_read_count_ch

    output:
    val(input_count) into input_check_ch

    script:
    """
    #!/usr/bin/env python
    for requested_count in $target_reads:
        if int(requested_count) >= $input_count:
            raise ValueError(f'Requested number of reads ({requested_count}) is greater than the input ({$input_count})')
    """
}


// TODO - run depthofcoverage too?
// TODO - total reads-duplicates - need to mark duplicates first

/* ------------------------------------------------------------------------ *
 *   Downsample BAM
 * ------------------------------------------------------------------------ */

// calculate percent to downsample to
process calc_percent_downsample {
    conda "$baseDir/env/downsample.yml"
    tag "$depth"

    input:
    val(depth) from depths_ch
    val(input_count) from input_check_ch

    output:
    tuple val(depth),
          file("${depth}_percent.txt") into downsample_percent_ch

    script:
    """
    #!/usr/bin/env python
    target = int('$depth'.strip('reads')) / $input_count
    with open('${depth}_percent.txt', 'w') as f:
        f.write(str(target))
    """
}


// downsample BAM
process downsample {
    conda "$baseDir/env/downsample.yml"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"
    tag "${depth}_${rep}"

    input:
    tuple val(depth),
          file(percent_file) from downsample_percent_ch
    each(rep) from reps_ch

    output:
    tuple val(depth),
          val(rep),
          file("${depth}_${rep}.bam"), 
          file("${depth}_${rep}.downsample_metrics") into downsampled_bams_ch
    file("${depth}_${rep}.bai")

    script:
    """
    percent=\$(cat $percent_file)
    picard DownsampleSam \
        I=$params.input_bam \
        O=${depth}_${rep}.bam \
        M=${depth}_${rep}.downsample_metrics \
        CREATE_INDEX=true \
        RANDOM_SEED=null \
        P=\$percent
    """
}


// calculate total reads with duplicates removed
process mark_dups {
    conda "$baseDir/env/downsample.yml"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"
    tag "${depth}_${rep}"

    input:
    tuple val(depth), 
          val(rep), 
          file(bam),
          file(bam_metrics) from downsampled_bams_ch

    output:
    file("${depth}_${rep}.rmdup.bam")
    file("${depth}_${rep}.rmdup.bai")
    tuple val(depth),
          val(rep),
          file("${depth}_${rep}.rmdup_metrics"),
          file(bam_metrics) into rmdup_metrics_ch

    script:
    """
    picard MarkDuplicates \
        I=$bam \
        O=${depth}_${rep}.rmdup.bam \
        M=${depth}_${rep}.rmdup_metrics \
        CREATE_INDEX=true
    """
}


// TODO - run depthofcoverage?


// combine metrics files into one per sample
process combine_metrics {
    tag "${depth}_${rep}"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"

    input:
    tuple val(depth),
          val(rep),
          file(rmdup),
          file(ds) from rmdup_metrics_ch

    output:
    file("${depth}_${rep}_metrics.csv") into combined_metrics_ch

    script:
    """
    # gather metics from downsampling metrics file
    ds_total_reads=\$(head -n8 $ds | tail -n1 | cut -f1)
    ds_pf_reads=\$(head -n8 $ds | tail -n1 | cut -f2)
    ds_read_length=\$(head -n8 $ds | tail -n1 | cut -f3)
    ds_total_bases=\$(head -n8 $ds | tail -n1 | cut -f4)
    ds_pf_bases=\$(head -n8 $ds | tail -n1 | cut -f5)
    
    # gather metrics from rmdup metics file
    rmdup_unpaired_reads_examined=\$(head -n8 $rmdup | tail -n1 | cut -f2)
    rmdup_read_pairs_examined=\$(head -n8 $rmdup | tail -n1 | cut -f3)
    rmdup_secondary_or_supplementary_reads=\$(head -n8 $rmdup | tail -n1 | cut -f4)
    rmdup_unmapped_reads=\$(head -n8 $rmdup | tail -n1 | cut -f5)
    rmdup_unpaired_read_duplicates=\$(head -n8 $rmdup | tail -n1 | cut -f6)
    rmdup_read_pair_duplicates=\$(head -n8 $rmdup | tail -n1 | cut -f7)
    rmdup_read_pair_optical_duplicates=\$(head -n8 $rmdup | tail -n1 | cut -f8)
    rmdup_percent_duplication=\$(head -n8 $rmdup | tail -n1 | cut -f9)
    rmdup_estimated_library_size=\$(head -n8 $rmdup | tail -n1 | cut -f10)

    # print metrics to file
    echo -e "${depth}_${rep},\\
    \$(echo ${depth} | tr -cd '0-9'),\\
    \$(echo ${rep} | tr -cd '0-9'),\\
    \$ds_total_reads,\\
    \$ds_pf_reads,\\
    \$ds_read_length,\\
    \$ds_total_bases,\\
    \$ds_pf_bases,\\
    \$rmdup_unpaired_reads_examined,\\
    \$rmdup_read_pairs_examined,\\
    \$rmdup_secondary_or_supplementary_reads,\\
    \$rmdup_unmapped_reads,\\
    \$rmdup_unpaired_read_duplicates,\\
    \$rmdup_read_pair_duplicates,\\
    \$rmdup_read_pair_optical_duplicates,\\
    \$rmdup_percent_duplication,\\
    \$rmdup_estimated_library_size" > ${depth}_${rep}_metrics.csv
    """
}


// TODO - add input BAM
// combine metrics files from all samples into one
process combine_samples {
    publishDir "${params.outdir}/downsample/", mode: "copy"

    input:
    file(metrics_file) from combined_metrics_ch.collect()

    output:
    file("combined_metrics.csv")

    script:
    """
    # make header
    echo -e "sample,\\
    target_reads,\\
    rep,\\
    ds_total_reads,\\
    ds_pf_reads,\\
    ds_read_length,\\
    ds_total_bases,\\
    ds_pf_bases,\\
    rmdup_unpaired_reads_examined,\\
    rmdup_read_pairs_examined,\\
    rmdup_secondary_or_supplementary_reads,\\
    rmdup_unmapped_reads,\\
    rmdup_unpaired_read_duplicates,\\
    rmdup_read_pair_duplicates,\\
    rmdup_read_pair_optical_duplicates,\\
    rmdup_percent_duplication,\\
    rmdup_estimated_library_size" > combined_metrics.csv

    # add sample
    cat $metrics_file >> combined_metrics.csv
    """
}
