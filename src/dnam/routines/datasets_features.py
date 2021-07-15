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


def get_status_dict(dataset: str):
    d = None
    if dataset == "GSE42861":
        d = {"Control": "Normal", "Case": "Patient"}
    elif dataset == "GSE80417":
        d = {"Control": 1, "Case": 2}
    elif dataset == "GSE84727":
        d = {"Control": 1, "Case": 2}
    elif dataset == "GSE125105":
        d = {"Control": "control", "Case": "case"}
    elif dataset == "GSE147221":
        d = {"Control": "Control", "Case": "Case"}
    elif dataset == "GSE152027":
        d = {"Control": "CON", "Case": "SCZ"}
    elif dataset == "GSE168739":
        d = {"Control": "", "": ""}
    elif dataset == "GSE111629":
        d = {"Control": "PD-free control", "Case": "Parkinson's disease (PD)"}
    elif dataset == "GSE128235":
        d = {"Control": "control", "Case": "case"}
    elif dataset == "GSE72774":
        d = {"Control": "control", "Case": "PD"}
    elif dataset == "GSE53740":
        d = {"Control": "Control", "Case": "FTD"}
    elif dataset == "GSE87648":
        d = {"Control": "HL", "Case": "CD"}
    elif dataset == "GSE144858":
        d = {"Control": "control", "Case": "Alzheimer's disease"}
    return d


def get_status_case_name(dataset: str):
    d = {
        'GSE42861': 'Rheumatoid Arthritis',
        'GSE80417': 'Schizophrenia',
        'GSE84727': 'Schizophrenia',
        'GSE125105': 'Depression',
        'GSE147221': 'Schizophrenia',
        'GSE152027': 'Schizophrenia',
        'GSE168739': '',
        'GSE111629': "Parkinson's disease",
        'GSE128235': 'Depression',
        'GSE72774': "Parkinson's disease",
        'GSE53740': 'Frontotemporal dementia',
        'GSE87648': "Crohn's disease",
        'GSE144858': "Alzheimer's disease",
    }
    return d[dataset]


def get_sex_dict(dataset: str):
    d = None
    if dataset == "GSE42861":
        d = {"F": "f", "M": "m"}
    elif dataset == "GSE80417":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE84727":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE125105":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE147221":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE152027":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE168739":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE111629":
        d = {"F": "Female", "M": "Male"}
    elif dataset == "GSE128235":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE72774":
        d = {"F": "female", "M": "male"}
    elif dataset == "GSE53740":
        d = {"F": "FEMALE", "M": "MALE"}
    elif dataset == "GSE87648":
        d = {"F": "F", "M": "M"}
    elif dataset == "GSE144858":
        d = {"F": "Female", "M": "Male"}
    return d
