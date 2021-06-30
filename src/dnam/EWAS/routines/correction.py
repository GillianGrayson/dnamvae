from statsmodels.stats.multitest import multipletests

def correct_pvalues(data: dict, metrics: list):

    for m in metrics:
        reject, pvals_corr, alphacSidak, alphacBonf = multipletests(data[m], 0.05, method='fdr_bh')
        data[f'{m}_fdr_bh'] = pvals_corr
        reject, pvals_corr, alphacSidak, alphacBonf = multipletests(data[m], 0.05, method='bonferroni')
        data[f'{m}_bonferroni'] = pvals_corr

    return data