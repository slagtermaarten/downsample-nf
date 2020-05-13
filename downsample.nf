// TODO - add check of input data
input = file("$params.input_bam")


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
  .fromList( params.target_reads )
  .map{ "${it}reads" }
  .set{ depths_ch }


// make channel for each repeat
Channel
  .from( 1 .. params.num_repeats )
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
    publishDir "${params.outdir}/$depth", mode: "copy"
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
