process macs3 {

    tag "${id}"
    label 'med_mem'

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( bamfile ), path( bamfile_ctrl, stageAs: "control/*" )

    output:
    tuple val( id ), path( "macs3/*_peaks.tsv" ), emit: excel
    tuple val( id ), path( "macs3/*.bed" ), emit: bed
    tuple val( id ), path( "macs3/*.bdg" ), emit: bg
    tuple val( id ), path( "macs3/*_cutoff_analysis.txt" ), emit: cutoffs

    script:
    """
    samtools merge -@ ${task.cpus} -o merged.bam ${bamfile}
    samtools view -bS merged.bam -h -F 20 -o f.bam
    samtools view -bS merged.bam -h -f 16 -o r.bam

    samtools view -bS "${bamfile_ctrl}" -h -F 20 -o f_ctrl.bam
    samtools view -bS "${bamfile_ctrl}" -h -f 16 -o r_ctrl.bam

    for f in {f,r}.bam
    do
        strand=\$(basename \$f .bam)
        ctrl=\$(basename \$f .bam)_ctrl.bam
        macs3 callpeak \
            --treatment "\$f" \
            --control "\$ctrl" \
            --format BAMPE \
            --call-summits \
            -p 0.2 \
            --cutoff-analysis \
            --outdir macs3 \
            --name "${id}.\$strand" \
            --bdg \
            --trackline \
            --gsize 4000000 \
            --buffer-size ${Math.round(task.memory.toMega() * 0.8 / 800)}
    done

    for f in macs3/*.xls
    do
        # these aren't actually xls files
        mv \$f \$(dirname \$f)/\$(basename \$f .xls).tsv
    done

    rm merged.bam

    """
}