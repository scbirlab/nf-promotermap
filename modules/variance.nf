process get_peak_variance {

    tag "${id}"

    publishDir( 
        "${params.outputs}/peaks", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( coverage )

    output:
    tuple val( id ), path( "peak-cov-rpkm-var.tsv" )

    script:
    """
    #!/usr/bin/env python

    import numpy as np
    import pandas as pd

    df = pd.read_csv("${coverage}", sep="\\t")
    df = df.set_index(df.columns.tolist()[:6])
    df.columns = df.columns.map(lambda x: x.split(".sorted.bam")[0])
    bam_cols = sorted(df.columns)
    df = df[bam_cols]
    print(df.head())

    peak_size_kb = pd.Series(
        np.abs(df.index.get_level_values("peak_start").values - df.index.get_level_values("peak_end").values) / 1000.,
        index=df.index,
    )
    print(peak_size_kb)
    bin_sums = df.sum(axis=0)
    df_rpkm = (df.values / peak_size_kb.values[:,None]) / bin_sums.values[None,:]
    df_rpkm = pd.DataFrame(df_rpkm, index=df.index, columns=df.columns)

    df2 = (
        df_rpkm
        .assign(
            peak_sum=df_rpkm.sum(axis=1),
            peak_var=df_rpkm.var(axis=1),
            peak_max=df_rpkm.idxmax(axis="columns"),
        )
    )
    df2.to_csv("peak-cov-rpkm-var.tsv", sep="\\t")

    """
}