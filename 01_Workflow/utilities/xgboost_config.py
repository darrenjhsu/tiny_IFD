
def xgboost_config():
    c = {}
    c['n_estimators'] = 100
    c['max_depth'] = 6
    c['gamma'] = 0.3
    c['suffix'] = f'{c["n_estimators"]}_{c["max_depth"]}_{c["gamma"]}'


    c['norm_cols'] = [ # Normal columns
        'Dihedral', 'Elec', 'vdW', 'Solvent', 
        'ligBond', 'ligAngle', 'ligDihe', 
        'ligvdW', 'ligvdW14', 'extravdW', #'extravdW14',
        'apovdW', 'apovdW14',
        # 'holovdW', 'holovdW14',
        'ligelec', 'ligelec14', 'extraelec', #'extraelec14',
        'apoelec', 'apoelec14',
        # 'holoelec', 'holoelec14',
        'changeSASA', 'changeCoreSASA', 
        'HBond', 'Contact',
        'proteinRMSD', 'proteinCoreRMSD', 
        'templateOverlap', 'temporalRMSD',
    ]

    c['norm_cols'] = [ 
    'Dihedral', 'Elec', 'vdW', 'Solvent', 
    'ligBond', 'ligAngle', 'ligDihe', 
    'ligvdW', 'ligvdW14', 'extravdW',
    'apovdW', 'apovdW14',
    'ligelec', 'ligelec14', 'extraelec',
    'apoelec', 'apoelec14',
    'allChangeSASA', 'ligandChangeSASA', 'ligandCSPSASA', 'ligandCHSASA', 'ligandPolarSASA',
    'HBond', 'Contact',
    'proteinRMSD', 'proteinCoreRMSD', 
    'templateOverlap', 'temporalRMSD',
]

#    c['norm_cols'] = [ # Simple columns
#        'vdW', 'ligvdW', 'extravdW',
#        'extraelec',
#        'changeSASA', 
#        'Contact',
#        'templateOverlap'
#    ]
    
    c['glob_cols'] = ['proteinRMSD', 'proteinCoreRMSD', 'templateOverlap', 'temporalRMSD']
    c['glob_cols']= []
    c['raw_cols'] = ['proteinRMSD', 'proteinCoreRMSD', 'templateOverlap', 'temporalRMSD']
    c['raw_cols'] = []

    return c
