
def plot_aggregate_RMSD_one_ligand(MDR, key=None, mode='MD'):
    from matplotlib import pyplot as plt
    import numpy as np
#     MDR.gatherMagnitude()

    if key is None:
        if len(MDR.Ligands) == 1:
            key = list(MDR.Ligands.keys())[0]
            print(key)
        else:
            raise("MDR contains more than one ligand. Specify the ligand with key=<ligand> or use plot_aggregate_RMSD_per_ligand")

    fig, axs = plt.subplots(figsize=(10,6))
    RMSD_aggregate = []
    if type(mode) == list:
        for jj in MDR.Ligands[key].Poses:
            for mm in mode:
                for kk in MDR.Ligands[key].Poses[jj].traj[mm]:
                    if kk.hasRMSD:
                        RMSD_aggregate.append(kk.RMSD)
    else:    
        for jj in MDR.Ligands[key].Poses:
            for kk in MDR.Ligands[key].Poses[jj].traj[mode]:
                if kk.hasRMSD:
                    RMSD_aggregate.append(kk.RMSD)

    RMSD_aggregate = np.array(RMSD_aggregate).flatten()
    axs.hist(RMSD_aggregate,bins=np.linspace(1,10,50))
    axs.grid(alpha=0.3)

    plt.xlabel('RMSD (Å)')
    plt.ylabel('Counts')


    
def plot_aggregate_RMSD_per_ligand(MDR, mode='MD'):
    from matplotlib import pyplot as plt
    import numpy as np
    MDR.gatherMagnitude()
    fig, axs = plt.subplots(figsize=(18,(MDR.numSystems+5)//5+1),nrows=(MDR.numSystems+5)//5,ncols=5,sharex='col',sharey='row')
    for ii, key in enumerate(MDR.Ligands):
    #     print(ii)
        RMSD_aggregate = []
        for jj in MDR.Ligands[key].Poses:
            if mode == 'MD':
                for kk in MDR.Ligands[key].Poses[jj].trajMD:
                    if kk.hasRMSD:
                        RMSD_aggregate.append(kk.RMSD)
                ymax = MDR.frameMD // 500
            if mode == 'QR':
                for kk in MDR.Ligands[key].Poses[jj].trajQR:
                    if kk.hasRMSD:
                        RMSD_aggregate.append(kk.RMSD)
                ymax = MDR.frameQR // 500
            if mode == 'EM':
                for kk in MDR.Ligands[key].Poses[jj].trajEM:
                    if kk.hasRMSD:
                        RMSD_aggregate.append(kk.RMSD)            
                ymax = MDR.frameEM // 500
        RMSD_aggregate = np.array(RMSD_aggregate).flatten()
        axs[ii//5, ii%5].hist(RMSD_aggregate,bins=np.linspace(1,10,50))
        axs[ii//5, ii%5].grid(alpha=0.3)
        axs[ii//5, ii%5].text(1, 0.8 * ymax, f'{key}')
        axs[ii//5, ii%5].set_ylim([0,ymax])
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.00, wspace=0.00)
    plt.subplots_adjust(left=0.06)
    plt.subplots_adjust(bottom=0.06)
    fig.text(0.00, 0.5, 'Counts', ha='center', rotation='vertical')
    fig.text(0.5, 0.00, 'RMSD (Å)', va='center')

    
def plot_cumulative_RMSD(MDR, mode='MD'):
    from matplotlib import pyplot as plt
    import numpy as np
    RMSD_all = []
    for ii in MDR.Ligands:
        for jj in MDR.Ligands[ii].Poses:
            if type(mode) == list:
                for mm in mode:
                    for kk in MDR.Ligands[ii].Poses[jj].traj[mm]:
                        if kk.hasRMSD:
                            RMSD_all.append(kk.RMSD)
            else:
                for kk in MDR.Ligands[ii].Poses[jj].traj[mode]:
                    if kk.hasRMSD:
                        RMSD_all.append(kk.RMSD)
    RMSD_all = np.array(RMSD_all).flatten()
    cumulative_threshold = np.linspace(0,10,200)
    cumulative_rmsd = np.zeros_like(cumulative_threshold)
    for idx, ii in enumerate(cumulative_threshold):
        cumulative_rmsd[idx] = np.sum((RMSD_all < ii)) / len(RMSD_all)
    plt.figure(figsize=(10,6))
    plt.plot(cumulative_threshold, cumulative_rmsd)
    plt.grid()
    plt.xlabel('RMSD (Å)',fontsize=20)
    plt.ylabel('Cumulative portion of frames',fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    return RMSD_all
    
    
def plot_histogram_RMSD(MDR, mode='MD', density=True):
    from matplotlib import pyplot as plt
    import numpy as np
    RMSD_all = []
    for ii in MDR.Ligands:
        for jj in MDR.Ligands[ii].Poses:
            if type(mode) == list:
                for mm in mode:
                    for kk in MDR.Ligands[ii].Poses[jj].traj[mm]:
                        if kk.hasRMSD:
                            RMSD_all.append(kk.RMSD)
            else:
                for kk in MDR.Ligands[ii].Poses[jj].traj[mode]:
                    if kk.hasRMSD:
                        RMSD_all.append(kk.RMSD)
    RMSD_all = np.array(RMSD_all).flatten()
    plt.figure(figsize=(10,6))
    try:
        plt.hist(RMSD_all,bins=np.linspace(1,10,MDR.frame['MD']//5000),density=density)
    except:
        plt.hist(RMSD_all,bins=np.linspace(1,10,271),density=density)
    plt.grid()
    plt.xlabel('RMSD (Å)',fontsize=20)
    plt.ylabel('Prob density of frames',fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    return RMSD_all

def gatherRMSD(MDR, mode='MD'):
    import numpy as np
    RMSD_all = []
    for ii in MDR.Ligands:
        for jj in MDR.Ligands[ii].Poses:
            if mode == 'init':
                RMSD_all.append(
                    MDR.Ligands[ii].Poses[jj].initialRMSD)
            if mode == 'MD':
                for kk in MDR.Ligands[ii].Poses[jj].trajMD:
                    if kk.hasRMSD:
                        RMSD_all.append(kk.RMSD)

            if mode == 'QR':
                for kk in MDR.Ligands[ii].Poses[jj].trajQR:
                    if kk.hasRMSD:
                        RMSD_all.append(kk.RMSD)

            if mode == 'EM':
                for kk in MDR.Ligands[ii].Poses[jj].trajEM:
                    if kk.hasRMSD:
                        RMSD_all.append(kk.RMSD)            

    RMSD_all = np.array(RMSD_all).flatten()
    return RMSD_all