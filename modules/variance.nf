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
    from scipy.stats import skew, kurtosis

    def norm_df(df):
        return df.values / df.values.sum(axis=1, keepdims=True)

    def entropy(df):
        normed = norm_df(df)
        return -np.nansum(xlogy(normed , normed), axis=1)


    def sparsity_index(df):
        n = df.shape[1]
        sqrt_n = np.sqrt(n)
        return (
            sqrt_n - np.linalg.norm(df.values, 1, axis=1, keepdims=True) 
            / np.linalg.norm(df.values, 2, axis=1, keepdims=True)
        ) / (sqrt_n - 1.)


    def _var(df):
        v = np.arange(df.shape[1], dtype=np.float64)[None]
        normed = norm_df(df)
        means = (normed * v).sum(axis=1, keepdims=True)
        means2 = (normed * np.square(v)).sum(axis=1, keepdims=True)
        return (means2 - np.square(means)).flatten()


    def _bc(df):
        normed = norm_df(df)
        g = skew(normed, axis=1, keepdims=True)
        k = kurtosis(normed, fisher=False, axis=1, keepdims=True) 
        n = normed.shape[1]
        return (np.square(g) + 1.) / (k + 3. * np.square(n - 1.)) / ((n - 2.) * (n - 3.))

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
            peak_var=_var(df_rpkm),
            peak_entropy=entropy(df_rpkm),
            peak_sparsity=sparsity_index(df_rpkm),
            peak_bc=_bc(df_rpkm),
            peak_max=df_rpkm.idxmax(axis="columns"),
            peak_max_n=lambda x: x["peak_max"].map(sorted(df_rpkm.columns.tolist()).index),
        )
    )
    df2.to_csv("peak-cov-rpkm-var.tsv", sep="\\t")

    """
}