process bam2wig {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/bigwig", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "coverage.bw" )

    script:
    """
    samtools index "${bamfile}"
    bamCoverage \
        -b "${bamfile}" \
        --outFileName coverage.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors ${task.cpus}
    
    """
}