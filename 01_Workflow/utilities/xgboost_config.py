
def xgboost_config():
    c = {}
    c['n_estimators'] = 100
    c['max_depth'] = 6
    c['gamma'] = 0.2
    c['additional_suffix'] = 'allcol_notempo_pharma'
    
    c['suffix'] = f'{c["n_estimators"]}_{c["max_depth"]}_{c["gamma"]}'
    
    c['suffix'] += c['additional_suffix']


    c['norm_cols'] = [ 
    'Dihedral', 'Elec', 'vdW', 'Solvent', 
    'ligBond', 'ligAngle', 'ligDihe', 
    'ligvdW', 'ligvdW14', 'extravdW',
    'apovdW', 'apovdW14',
    'ligelec', 'ligelec14', 'extraelec',
    'apoelec', 'apoelec14',
    'allChangeSASA', 'ligandChangeSASA', 'ligandCSPSASA', 'ligandCHSASA', 'ligandPolarSASA',
    'HBond', 'Contact', # Contact is the same as C-C below
    'proteinRMSD', 'proteinCoreRMSD', 
    'templateOverlap', 'pharmacoOverlap',#'temporalRMSD',
    'H-H', 'H-C', 'H-N',
    'H-O', 'H-S', 'H-P', 'H-halo', 'C-H', 'C-C', 'C-N', 'C-O', 'C-S', 'C-P',
    'C-halo', 'N-H', 'N-C', 'N-N', 'N-O', 'N-S', 'N-P', 'N-halo', 'O-H',
    'O-C', 'O-N', 'O-O', 'O-S', 'O-P', 'O-halo', 'S-H', 'S-C', 'S-N', 'S-O',
    'S-S', 'S-P', 'S-halo'
] # This is the original version which has a number of useless columns such as 'C-P' etc.

    c['norm_cols'] = [ 
    'vdW', 'Solvent',  #From output
    'ligvdW', 'extravdW', #From cpptraj
    'ligelec', 'ligelec14', 'extraelec', #Derived from rom cpptraj
    'allChangeSASA', 'ligandChangeSASA', 'ligandCSPSASA', 'ligandCHSASA', 'ligandPolarSASA', # from SASA
    'HBond', # from getHBHC
    'proteinRMSD', 'proteinCoreRMSD', # from MDR processing
    'templateOverlap', 'pharmacoOverlap', # from MDR processing
    'H-H', 'H-C', 'H-N', 'H-O', 'H-S', 'H-halo',
    'C-H', 'C-C', 'C-N', 'C-O', 'C-S', 'C-halo', 
    'N-H', 'N-C', 'N-N', 'N-O', 'N-S', 'N-halo', 
    'O-H', 'O-C', 'O-N', 'O-O', 'O-S', 'O-halo', 
    'S-H', 'S-C', 'S-N', 'S-O', 'S-halo' # Contacts are from MDR processing
]
    c['glob_cols']= []
    c['raw_cols'] = []


    return c


def create_config_list(csv, limit=None):

    import pandas as pd

    df = pd.read_csv(csv)
    print(f'We will have {len(df)} conditions to try')

    if limit is not None:
        df = df[:limit]

    df_list = df.to_dict(orient='records')

    c = {}
    c['norm_cols'] = xgboost_config()['norm_cols']
    c['glob_cols']= []
    c['raw_cols'] = []
    
    c_list = []
    for r in df_list:
        for key in r.keys():
            c[key] = r[key]
        c['suffix'] = f'{c["n_estimators"]}_{c["max_depth"]}_{c["gamma"]}_{c["colsample_bytree"]}_{c["subsample"]}'
        c_list.append(c.copy())

    return c_list
        

