import pickle
import numpy as np
import glob, os
import pandas as pd
import os, sys, gc, glob



def normalize(df, return_reference=False):
    df_mean = df.mean(axis=0)
    df_std = df.std(axis=0)
    if not return_reference:
        return (df-df_mean) / df_std
    else:
        return (df-df_mean) / df_std, df_mean, df_std
def normalize_to_ref(df, ref_mean, ref_std):
    return (df-ref_mean) / ref_std
def concentrate(df, rmsd, thres, ratio=0.1):#, do_nothing_if_over = True):
    assert len(df) == len(rmsd), "df must be of same length to rmsd"
    rmsd_sel = rmsd <= thres
    if np.sum(rmsd_sel) / len(rmsd_sel) > ratio:
        return df, rmsd
    else:
        df2 = df
        df2['rmsd'] = rmsd
        good = df2.loc[df['rmsd'] <= thres]
        bad = df2.loc[~(df['rmsd'] < thres)].sample(int(1/ratio-1)*len(good))
        df3 = pd.concat([good, bad])
        rmsd3 = df3['rmsd']
        df3.drop(['rmsd'], axis=1, inplace=True)
        return df3, rmsd3


def clean_train_data(X_raw, y_raw = None, std_threshold = 5, rmsd_threshold = 15, dataset='train'):
    from xgboost_config import xgboost_config
    norm_cols = xgboost_config()['norm_cols']
    glob_cols = xgboost_config()['glob_cols']
    train_mask = ~np.any(np.abs(X_raw[norm_cols]) > std_threshold, axis=1) & ~np.any(np.abs(X_raw[glob_cols]) > std_threshold, axis=1)
    if dataset == 'train':
        train_mask = train_mask & (y_raw < rmsd_threshold)
        yv = y_raw[train_mask]
        y = yv < 2.5
    elif dataset == 'val':
        yv = y_raw[train_mask]
        y = yv < 2.5
    elif dataset == 'test':
        yv = y = None
    else:
        raise NotImplementedError(f'dataset must be one of train, val, or test, not "{dataset}"')
    X = X_raw[train_mask]

    return X, yv, y, train_mask 

