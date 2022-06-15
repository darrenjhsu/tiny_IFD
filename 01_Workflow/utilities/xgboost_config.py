
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
    'HBond', 'Contact',
    'proteinRMSD', 'proteinCoreRMSD', 
    'templateOverlap', 'pharmacoOverlap',#'temporalRMSD',
    'H-H', 'H-C', 'H-N',
    'H-O', 'H-S', 'H-P', 'H-halo', 'C-H', 'C-C', 'C-N', 'C-O', 'C-S', 'C-P',
    'C-halo', 'N-H', 'N-C', 'N-N', 'N-O', 'N-S', 'N-P', 'N-halo', 'O-H',
    'O-C', 'O-N', 'O-O', 'O-S', 'O-P', 'O-halo', 'S-H', 'S-C', 'S-N', 'S-O',
    'S-S', 'S-P', 'S-halo'
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

    print(c)

    return c
