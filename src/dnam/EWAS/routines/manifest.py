import pandas as pd
import os


def process_str_elem(x, delimiter: str = ';'):
    if isinstance(x, str):
        elems = x.split(';')
        elems = list(set(elems))
        elems = delimiter.join(elems)
    else:
        elems = ''
    return elems


def get_manifest(platform="GPL13534"):

    fn_pkl = f"E:/YandexDisk/Work/pydnameth/datasets/{platform}/manifest/manifest.pkl"
    if os.path.isfile(fn_pkl):
        manifest = pd.read_pickle(fn_pkl)
    else:
        fn = f"E:/YandexDisk/Work/pydnameth/datasets/{platform}/manifest/manifest.xlsx"
        manifest = pd.read_excel(fn, index_col="CpG")
        manifest['Gene'] = manifest['Gene'].apply(process_str_elem)
        manifest.to_pickle(fn_pkl)

    return manifest
