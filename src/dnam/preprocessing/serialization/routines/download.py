import pandas as pd
import GEOparse
from tqdm import tqdm

def download_betas_and_pvals_from_gsms(gsms, path):

    for gsm_id, gsm in tqdm(enumerate(gsms)):
        while True:
            try:
                gsm_data = GEOparse.get_GEO(geo=gsm, destdir= f"{path}/gsms", include_data=True, how="full", silent=True)
            except ValueError:
                continue
            except ConnectionError:
                continue
            except IOError:
                continue
            break

        gsm_data.table.set_index('ID_REF', inplace=True)
        if gsm_id == 0:
            betas = pd.DataFrame(0, index=gsm_data.table.index, columns=gsms)
            pvals = pd.DataFrame(0, index=gsm_data.table.index, columns=gsms)
        betas[gsm] = gsm_data.table['VALUE']
        pvals[gsm] = gsm_data.table['DetectionPval']

    betas = betas.T
    betas.index.name = "subject_id"
    pvals = pvals.T
    betas.index.name = "subject_id"

    print(f"Number of NaNs in betas: {betas.isna().values.sum()}")
    print(f"Number of NaNs in pvals: {pvals.isna().values.sum()}")

    return betas, pvals
