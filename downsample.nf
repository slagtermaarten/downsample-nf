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
input = file("$params.input_bam", checkIfExists: true)

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
    file(input_bam) from input

    output:
    file("input_bam_read_count.txt") into input_bam_read_count_ch

    script:
    """
    samtools view -c -F 260 $input_bam > input_bam_read_count.txt
    """
}


// are there enough reads to meet num_reads? - if yes, input to downsample process 
// TODO? - check inputs are integers
process input_bam_qc {
    conda "env/downsample.yml"
    publishDir "${params.outdir}", mode: "copy"

    input:
    file(input_count) from input_bam_read_count_ch

    output:
    val(true) into input_check_ch

    script:
    """
    #!/usr/bin/env python
    # get input read count
    with open("$input_count") as f:
        read_count = int(f.readline())
    
    # compare to num_reads, exit if requested num is larger than actual number
    for requested_count in $params.num_reads:
        if int(requested_count) >= read_count:
            raise ValueError(f'Requested number of reads ({requested_count}) is greater than the input ({read_count})')
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
  .set{combined_ch_1}


// TODO - calculate percent to downsample to

// downsample BAM
process downsample {
    conda "env/downsample.yml"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"
    tag "${depth}_${rep}"

    input:
    val(flag) from input_check_ch
    set(depth, rep) from combined_ch_1

    output:
    file("${depth}_${rep}.bam") into downsampled_bams_ch
    file("${depth}_${rep}.bai")
    file("${depth}_${rep}.downsample_metrics")
    set(depth, rep) into combined_ch_2

    script:
    """
    picard DownsampleSam \
      I=$params.input_bam \
      O=${depth}_${rep}.bam \
      P=0.2 \
      M=${depth}_${rep}.downsample_metrics \
      RANDOM_SEED=null \
      CREATE_INDEX=true
    """
}


// TODO - check num of reads afterwards to make sure it worked okay??
//   - this duplcates metrics file from above - not needed?
/*
process output_qc {
  conda "env/downsample.yml"
  publishDir "${params.outdir}/downsample/$depth", mode: "copy"
  tag "${depth}_${rep}"

  input:
  file(bam) from downsampled_bams_ch
  set depth, rep from combined_ch_2

  output:
  file("${depth}_${rep}.read_count")

  script:
  """
  samtools view -c $bam > ${depth}_${rep}.read_count
  """
}
*/


// calculate total reads with duplicates removed
process mark_dups {
  conda "env/downsample.yml"
  publishDir "${params.outdir}/downsample/$depth", mode: "copy"
  tag "${depth}_${rep}"

  input:
  file(bam) from downsampled_bams_ch
  set(depth, rep) from combined_ch_2

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
