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
    tuple val( id ), path( "peaks.bed" ), emit: peaks
    tuple val( id ), path( "summits.bed" ), emit: summits
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

    samtools view -H "${bamfile_ctrl}" \
    | awk '\
        \$1 == "@SQ" { 
            for(i=1; i<=NF; i++) { 
                if(\$i ~ /^LN:/){
                    split(\$i, a, ":"); 
                    len+=a[2]
                }
            } 
        }
        END { print len }
    ' \
    > genome-size.txt

    for f in {f,r}.bam
    do
        strand=\$(basename \$f .bam)
        ctrl=\$(basename \$f .bam)_ctrl.bam
        macs3 callpeak \
            --treatment "\$f" \
            --control "\$ctrl" \
            --format BAMPE \
            --call-summits \
            -p 0.1 \
            --cutoff-analysis \
            --outdir macs3 \
            --name "${id}.\$strand" \
            --bdg \
            --trackline \
            --gsize \$(cat genome-size.txt) \
            --buffer-size ${Math.round(task.memory.toMega() * 0.8 / 800)}
    done

    for f in macs3/*.xls
    do
        # these aren't actually xls files
        mv \$f \$(dirname \$f)/\$(basename \$f .xls).tsv
    done

    summit_BED=(macs3/*.bed)
    head -n1 "\${summit_BED[0]}" \
    | sed 's/${id}\\.f/${id}/g' \
    | cat - \
        <(
            cat \
                <(awk -v OFS='\\t' 'NR>1 { print \$0, "+" }' "\${summit_BED[0]}") \
                <(awk -v OFS='\\t' 'NR>1 { print \$0, "-" }' "\${summit_BED[1]}") \
            | sort -k1,1 -k2,2n
        ) \
    > "summits.bed"

    peak_tsv=(macs3/*_peaks.tsv)
    cat \
        <(echo 'track name="${id} (peaks)" description="Peaks for ${id} (Made with MACS v3, 09/03/25)" visibility=1') \
        <(
            cat "\${peak_tsv[@]}" \
            | grep -v -e '^#' -e '^\$' -e '^chr\\s' \
            | awk -v OFS='\\t' '
                {
                    chr=\$1; start=\$2; end=\$3; name=\$10;
                    # infer strand from name tag
                    strand=".";
                    if (name ~ /\\.f_/) strand="+";
                    else if (name ~ /\\.r_/) strand="-";
                    # convert start to 0-based; end stays as-is (end-exclusive in BED)
                    print chr, start-1, end, name, 0, strand
                }
            ' \
            | sort -k1,1 -k2,2n \
        ) \
    > "peaks.bed"

    rm merged.bam

    """
}