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
    from scipy.special import xlogy

    def entropy(df):
        normed = df_rpkm.values / df.values.sum(axis=1, keepdims=True)
        return -np.nansum(xlogy(normed , normed), axis=1)

    df = pd.read_csv("${coverage}", sep="\\t")
    df = df.set_index(df.columns.tolist()[:6])
    df.columns = df.columns.map(lambda x: x.split(".sorted.bam")[0])
    bam_cols = sorted(df.columns)
    df = df[bam_cols]
    print(df.head())

    peak_size_kb = np.abs(df.index.get_level_values("peak_start").values - df.index.get_level_values("peak_end").values) / 1000.
    df_rpkm = df.values / (peak_size_kb[:,None] * df.sum(axis=0).values[None,:] / 1_000_000.)
    df_rpkm = pd.DataFrame(df_rpkm, index=df.index, columns=df.columns)

    df2 = (
        df_rpkm
        .assign(
            peak_sum=df_rpkm.sum(axis=1),
            peak_var=df_rpkm.var(axis=1),
            peak_entropy=entropy(df_rpkm),
            peak_max=df_rpkm.idxmax(axis="columns"),
        )
    )
    df2.to_csv("peak-cov-rpkm-var.tsv", sep="\\t")

    """
}