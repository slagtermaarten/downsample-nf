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

// check for required input data
if ( !params.input_bam )   { log.error "ERROR: input_bam is required"; exit 1 }
if ( !params.num_repeats ) { log.error "ERROR: num_repeats is required"; exit 1 }
if ( !params.num_reads )   { log.error "ERROR: num_reads is required"; exit 1 }

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
summary["Reads"]       = params.num_reads
summary["Num repeats"] = params.num_repeats
summary["Output dir"]  = params.outdir

log.info """
==============================================================================
Inputs:
${ summary.collect { k, v -> "  ${k.padRight(12)}: $v"}.join("\n") }
==============================================================================
"""


/* ------------------------------------------------------------------------ *
 *   Pipeline - check input BAM
 * ------------------------------------------------------------------------ */

// TODO - limit to bed region


// count total reads in input
process get_input_read_count {
    conda "env/downsample.yml"

    input:
    file(bam) from input_bam

    output:
    env(READ_COUNT) into input_bam_read_count_ch

    script:
    """
    READ_COUNT=\$(samtools view -c -F 260 $bam)
    """
}


// are there enough reads to meet num_reads? - if yes, input to downsample process 
// TODO? - check inputs are integers
process input_bam_qc {
    conda "env/downsample.yml"
    publishDir "${params.outdir}", mode: "copy"

    input:
    val(input_count) from input_bam_read_count_ch

    output:
    val(input_count) into input_check_ch

    script:
    """
    #!/usr/bin/env python
    for requested_count in $params.num_reads:
        if int(requested_count) >= $input_count:
            raise ValueError(f'Requested number of reads ({requested_count}) is greater than the input ({$input_count})')
    """
}


// TODO - run depthofcoverage too?
// TODO - total reads-duplicates - need to mark duplicates first


/* ------------------------------------------------------------------------ *
 *   Pipeline - downsample BAM
 * ------------------------------------------------------------------------ */

// make channel for each target depth
Channel
    .from( params.num_reads.split(',') )
    .map{ "${it}reads" }
    .set{ depths_ch }


// make channel for each repeat
Channel
    .from( 1..params.num_repeats.toInteger() )
    .map{ "rep$it" }
    .set{ reps_ch }


// combine channels
depths_ch
    .combine(reps_ch)
    .set{combined_ch}


// calculate percent to downsample to
process calc_percent_downsample {
    conda "env/downsample.yml"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"

    input:
    tuple val(depth), val(rep) from combined_ch
    val(input_count) from input_check_ch

    output:
    tuple val(depth), val(rep), file("${depth}_percent.txt") into downsample_percent_ch

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
    conda "env/downsample.yml"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"
    tag "${depth}_${rep}"

    input:
    tuple val(depth), val(rep), file(percent_file) from downsample_percent_ch

    output:
    tuple val(depth), val(rep), file("${depth}_${rep}.bam") into downsampled_bams_ch
    file("${depth}_${rep}.bai")
    file("${depth}_${rep}.downsample_metrics")

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
    conda "env/downsample.yml"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"
    tag "${depth}_${rep}"

    input:
    tuple val(depth), val(rep), file(bam) from downsampled_bams_ch

    output:
    file("${depth}_${rep}.rmdup.bam")
    file("${depth}_${rep}.rmdup.bai")
    file("${depth}_${rep}.rmdup_metrics")

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
