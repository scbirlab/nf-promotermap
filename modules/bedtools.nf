process gff2bed {

    tag "${id}"

    publishDir( 
        "${params.outputs}/genome", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( gff )

    output:
    tuple val( id ), path( "genes.bed" )

    script:
    """
    grep -v '^#' "${gff}" \
    | awk -v OFS='\\t' '
        \$3 == "gene" {
            chr=\$1; start=\$4-1; end=\$5; strand=\$7; attr=\$9;
            name="";
            n=split(attr, a, ";");
            for(i=1; i<=n; i++) {
                if(a[i] ~ /^Name=/) { sub(/^Name=/, "", a[i]); name=a[i] }
                else if(a[i] ~ /^old_locus_tag=/) { sub(/^old_locus_tag=/, "", a[i]); if(name == "") name=a[i]}
                else if(a[i] ~ /^locus_tag=/) { sub(/^locus_tag=/, "", a[i]); if(name == "") name=a[i]}
                else if(a[i] ~ /^gene=/) { sub(/^gene=/, "", a[i]); if(name == "") name=a[i]}
                else if(a[i] ~ /^ID=/) { sub(/^ID=/, "", a[i]); if(name == "") name=a[i] }
            }
            print chr, start, end, name, 0, strand
        }
        ' \
    | sort -k1,1 -k2,2n \
    > genes.bed

    """
}

process annotate_nearest_gene {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( summits ), path( ann_bed )

    output:
    tuple val( id ), path( "summits-ann.tsv" )

    script:
    """
    grep -v '^track' "${summits}" \
    | sort -k1,1 -k2,2n \
    > summits.sorted.bed

    printf 'chr\\tsummit_start\\tsummit_end\\tpeak_name\\tscore\\tstrand\\tann_chr\\tann_start\\tann_end\\tann_name\\tann_score\\tann_strand\\tann_offset\\n' \
    > summits-ann.tsv
    bedtools closest -a summits.sorted.bed -b "${ann_bed}" \
        -s -D a -t first \
    >> summits-ann.tsv

    """
}


process extract_peak_sequences {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( peaks ), path( genome )

    output:
    tuple val( id ), path( "peak-seqs.fasta" ), emit: fasta
    tuple val( id ), path( "peak-seqs.tsv" ), emit: tsv

    script:
    """
    bedtools getfasta -s -fi "${genome}" -bed "${peaks}" -name > peak-seqs.fasta
    printf 'peak_id\\tpeak_name\\tchr\\tpeak_start\\tpeak_end\\tpeak_sequence\\n' \
    | cat - <(
        bedtools getfasta -s -fi "${genome}" -bed "${peaks}" -name -tab \
        | awk -F'\\t' -v OFS='\\t' '{ split(\$1, a, ":"); split(a[4], c0, "-"); split(c0[2], c1, "("); print \$1, a[1], a[3], c0[1], c1[1], \$2 }' \
    ) \
    > peak-seqs.tsv

    """
}

process get_peak_coverage {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( peaks ), path( bams )

    output:
    tuple val( id ), path( "peak-cov.bed" ), emit: bed
    tuple val( id ), path( "peak-cov.tsv" ), emit: tsv

    script:
    """
    for f in ${bams}
    do 
        samtools index -@ ${task.cpus} "\$f"
    done

    bedtools multicov -p -bams ${bams} -bed ${peaks} > peak-cov.bed
    printf 'chr\\tpeak_start\\tpeak_end\\tpeak_name\\tscore\\tstrand\\t${bams.join('\\t')}\\n' > peak-cov.tsv
    cat peak-cov.bed >> peak-cov.tsv

    """
}


process make_peak_summary_table {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( summits ), path( peaks )

    output:
    tuple val( id ), path( "peak-summary.tsv" )

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    (
        pd.read_csv("${summits}", sep="\\t")
        .merge(
            pd.read_csv("${peaks}", sep="\\t"),
            how="left",
        )
        .to_csv("peak-summary.tsv", sep="\\t", index=False)
    )
    """
}
