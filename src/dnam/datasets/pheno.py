def get_status_pair(dataset: str):
    status = None
    if dataset == "GSE42861":
        status = ('subject', 'Status')
    elif dataset == "GSE80417":
        status = ('disease status', 'Status')
    elif dataset == "GSE84727":
        status = ('disease_status', 'Status')
    elif dataset == "GSE125105":
        status = ('diagnosis', 'Status')
    elif dataset == "GSE147221":
        status = ('status', 'Status')
    elif dataset == "GSE152027":
        status = ('status', 'Status')
    return status


def get_age_pair(dataset: str):
    age = None
    if dataset == "GSE42861":
        age = ("age", "Age")
    elif dataset == "GSE80417":
        age = ("age", "Age")
    elif dataset == "GSE84727":
        age = ("age", "Age")
    elif dataset == "GSE125105":
        age = ("age", "Age")
    elif dataset == "GSE147221":
        age = ("age", "Age")
    elif dataset == "GSE152027":
        age = ("ageatbloodcollection", "Age")
    return age


def get_sex_pair(dataset: str):
    sex = None
    if dataset == "GSE42861":
        sex = ("gender", "Sex")
    elif dataset == "GSE80417":
        sex = ("Sex", "Sex")
    elif dataset == "GSE84727":
        sex = ("Sex", "Sex")
    elif dataset == "GSE125105":
        sex = ("Sex", "Sex")
    elif dataset == "GSE147221":
        sex = ("Sex", "Sex")
    elif dataset == "GSE152027":
        sex = ("gender", "Sex")
    return sex

def get_status_vals_pairs(dataset: str):
    status_group = None
    if dataset == "GSE42861":
        status_group = [("Normal", "Status: Normal"), ("Patient", "Status: Rheumatoid Arthritis")]
    elif dataset == "GSE80417":
        status_group = [(1, "Status: Normal"), (2, "Status: Schizophrenia")]
    elif dataset == "GSE84727":
        status_group = [(1, "Status: Normal"), (2, "Status: Schizophrenia")]
    elif dataset == "GSE125105":
        status_group = [("control", "Status: Normal"), ("case", "Status: Depression")]
    elif dataset == "GSE147221":
        status_group = [("Control", "Status: Normal"), ("Case", "Status: Schizophrenia")]
    elif dataset == "GSE152027":
        status_group = [("CON", "Status: Normal"), ("SCZ", "Status: Schizophrenia")]
    return status_group


def get_sex_vals_pairs(dataset: str):
    status_group = None
    if dataset == "GSE42861":
        status_group = [("f", "F"), ("m", "M")]
    elif dataset == "GSE80417":
        status_group = [("F", "F"), ("M", "M")]
    elif dataset == "GSE84727":
        status_group = [("F", "F"), ("M", "M")]
    elif dataset == "GSE125105":
        status_group = [("F", "F"), ("M", "M")]
    elif dataset == "GSE147221":
        status_group = [("F", "F"), ("M", "M")]
    elif dataset == "GSE152027":
        status_group = [("F", "F"), ("M", "M")]
    return status_group
