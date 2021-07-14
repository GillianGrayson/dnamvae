def get_columns_dict(dataset: str):
    d = None
    if dataset == "GSE42861":
        d = {
            'Status': 'subject',
            'Age': 'age',
            'Sex': 'gender',
        }
    elif dataset == "GSE80417":
        d = {
            'Status': 'disease status',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE84727":
        d = {
            'Status': 'disease_status',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE125105":
        d = {
            'Status': 'diagnosis',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE147221":
        d = {
            'Status': 'status',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE152027":
        d = {
            'Status': 'status',
            'Age': 'ageatbloodcollection',
            'Sex': 'gender',
        }
    elif dataset == "GSE168739":
        d = {
            '': '',
            'Age': 'age',
            'Sex': 'gender',
        }
    elif dataset == "GSE111629":
        d = {
            'Status': 'disease state',
            'Age': 'age',
            'Sex': 'gender',
        }
    elif dataset == "GSE128235":
        d = {
            'Status': 'diagnosis',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE72774":
        d = {
            'Status': 'diseasestatus',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE53740":
        d = {
            'Status': 'diagnosis',
            'Age': 'age',
            'Sex': 'gender',
        }
    elif dataset == "GSE87648":
        d = {
            'Status': 'simplified_diagnosis',
            'Age': 'age',
            'Sex': 'Sex',
        }
    elif dataset == "GSE144858":
        d = {
            'Status': 'disease state',
            'Age': 'age',
            'Sex': 'Sex',
        }
    return d

def get_column_name(dataset: str, feature: str):
    d = get_columns_dict(dataset)
    return d[feature]

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
    elif dataset == "GSE168739":
        status = ('', '')
    elif dataset == "GSE111629":
        status = ('disease state', 'Status')
    elif dataset == "GSE128235":
        status = ('diagnosis', 'Status')
    elif dataset == "GSE72774":
        status = ('diseasestatus', 'Status')
    elif dataset == "GSE53740":
        status = ('diagnosis', 'Status')
    elif dataset == "GSE87648":
        status = ('simplified_diagnosis', 'Status')
    elif dataset == "GSE144858":
        status = ('disease state', 'Status')
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
    elif dataset == "GSE168739":
        age = ("age", "Age")
    elif dataset == "GSE111629":
        age = ("age", "Age")
    elif dataset == "GSE128235":
        age = ("age", "Age")
    elif dataset == "GSE72774":
        age = ("age", "Age")
    elif dataset == "GSE53740":
        age = ("age", "Age")
    elif dataset == "GSE87648":
        age = ("age", "Age")
    elif dataset == "GSE144858":
        age = ("age", "Age")
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
    elif dataset == "GSE168739":
        sex = ("gender", "Sex")
    elif dataset == "GSE111629":
        sex = ("gender", "Sex")
    elif dataset == "GSE128235":
        sex = ("Sex", "Sex")
    elif dataset == "GSE72774":
        sex = ("Sex", "Sex")
    elif dataset == "GSE53740":
        sex = ("gender", "Sex")
    elif dataset == "GSE87648":
        sex = ("Sex", "Sex")
    elif dataset == "GSE144858":
        sex = ("Sex", "Sex")
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
    elif dataset == "GSE168739":
        status_group = [("", ""), ("", "")]
    elif dataset == "GSE111629":
        status_group = [("PD-free control", "Status: Normal"), ("Parkinson's disease (PD)", "Status: Parkinson")]
    elif dataset == "GSE128235":
        status_group = [("control", "Status: Normal"), ("case", "Status: Depression")]
    elif dataset == "GSE72774":
        status_group = [("control", "Status: Normal"), ("PD", "Status: Parkinson")]
    elif dataset == "GSE53740":
        status_group = [("Control", "Status: Normal"), ("FTD", "Status: Frontotemporal dementia")]
    elif dataset == "GSE87648":
        status_group = [("HL", "Status: Normal"), ("CD", "Status: Crohn's disease")]
    elif dataset == "GSE144858":
        status_group = [("control", "Status: Normal"), ("Alzheimer's disease", "Alzheimer's disease")]
    return status_group


def get_sex_vals_pairs(dataset: str):
    status_group = None
    if dataset == "GSE42861":
        status_group = {"F": "f", "M": "m"}
    elif dataset == "GSE80417":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE84727":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE125105":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE147221":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE152027":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE168739":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE111629":
        status_group = {"F": "Female", "M": "Male"}
    elif dataset == "GSE128235":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE72774":
        status_group = {"F": "female", "M": "male"}
    elif dataset == "GSE53740":
        status_group = {"F": "FEMALE", "M": "MALE"}
    elif dataset == "GSE87648":
        status_group = {"F": "F", "M": "M"}
    elif dataset == "GSE144858":
        status_group = {"F": "Female", "M": "Male"}
    return status_group
