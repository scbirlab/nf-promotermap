process bam2wig {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/bigwig", 
        mode: 'copy',
    )

    input:
    tuple val( id ), path( bamfile )

    output:
    tuple val( id ), path( "${id}.bw" )

    script:
    """
    samtools index "${bamfile}"
    bamCoverage \
        -b "${bamfile}" \
        --outFileName "${id}.bw" \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors ${task.cpus}
    
    """
}


process plot_peaks {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/coverage", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bigwigs ), path( ann_bed )

    output:
    tuple val( id ), path( "matrix.mat.gz" ), emit: matrix
    tuple val( id ), path( "peak-heatmap.png" ), emit: plot

    script:
    """
    computeMatrix reference-point \
        --verbose \
        -S ${bigwigs} \
        --regionsFileName "${ann_bed}" \
        --referencePoint TSS \
        -a 500 \
        -b 500 \
        --binSize 10 \
        --skipZeros \
        --numberOfProcessors ${task.cpus} \
        -o matrix.mat.gz

    plotHeatmap \
        --verbose \
        --matrixFile matrix.mat.gz \
        --dpi 300 \
        --colorMap PRGn \
        --refPointLabel "Gene start" \
        --outFileName peak-heatmap.png
    
    """

}
