import pandas as pd


def get_forbidden_cpgs(path, types):
    forbidden_cpgs = set()
    for type in types:
        df = pd.read_excel(f"{path}/{type}.xlsx", header=None)
        tmp = df[0].values.tolist()
        forbidden_cpgs |= set(tmp)
    return list(forbidden_cpgs)
