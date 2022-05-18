import numpy as np

# Calculate contact for ODDT RF2, CPU version

ligand_map = {'A': 0, 'C': 0, 'N': 1, 'NA': 1, 'O': 2, 'OA': 2, 'F': 3, 'P': 4, 'S': 5, 'SA': 5, 'CL': 6,
              'BR': 7, 'I': 8}
protein_map = {'A': 0, 'C': 0, 'N': 1, 'NA': 1, 'O': 2, 'OA': 2, 'S': 3, 'SA': 3}


def contact(comp, cutoff=12.0, binsize=2.0, nbins=6):
    #comp is an mdtraj comp
    
    
    pro = comp.top.select("not residue 56 and not symbol H")
    proteinCoordinate = comp.xyz[0][pro]*10
    proteinType = np.array([protein_map[x.element.symbol.upper()] for x in np.array(list(comp.top.atoms))[pro]])
    lig = comp.top.select("residue 56 and not symbol H")
    ligandCoordinate = comp.xyz[0][lig]*10
    ligandType = np.array([ligand_map[x.element.symbol.upper()] for x in np.array(list(comp.top.atoms))[lig]])

    nfeatures = 216

    # Feature is ordered by ligtype, proteintype, bins
    features = np.zeros((1,nfeatures),dtype=int)
    
    # Calculate contact
    for ii in range(len(ligandType)):
        lc = ligandCoordinate[ii]
        lt = ligandType[ii]
        for jj in range(len(proteinType)):
            pc = proteinCoordinate[jj]
            pt = proteinType[jj]
            dist = np.sqrt(((pc - lc)**2).sum())
            if dist < 12.0:
#                 print(lc, lt, pc, pt, int(lt * 54 + pt * 6 + dist // 2))
#                 features[0,int(lt * 24 + pt * 6 + dist // 2)] += 1
                features[0,int(lt * 24 + pt * 6 + dist // 2)] += 1
    

    return features