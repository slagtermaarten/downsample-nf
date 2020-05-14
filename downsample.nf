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
input = file("$params.input_bam") // , checkIfExists: true  - TODO add this when using real data

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
 *   Pipeline
 * ------------------------------------------------------------------------ */

process input_bam {
    /*
     * TODO - check num reads in input bam and make sure there is enough - input to downsample process
     */
    publishDir "${params.outdir}", mode: "copy"
    input:
        file(input_bam) from input
    output:
        val(true) into input_check_ch
        file("summary.txt")
    script:
    """
    # check num reads in input BAM

    # add to summary text file
    echo input BAM XXXXX reads > summary.txt

    # check that there are enough reads to match target inputs
    """
}


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


// downsample BAM
process downsample {
    /*
     *  Downsample TODO - add actual downsampling script
     */
    tag "${depth}_${rep}"
    publishDir "${params.outdir}/downsample/$depth", mode: "copy"
    input:
        val flag from input_check_ch
        set depth, rep from combined_ch
    output:
        file "${depth}_${rep}.bam"
    script:
    """
    echo ${depth}_${rep} > ${depth}_${rep}.bam
    
    """
}

// TODO - index BAM??? - can add as flag to picard?

// TODO - check num of reads afterwards to make sure it worked okay??
