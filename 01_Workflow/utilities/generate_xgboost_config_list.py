from itertools import product
n_est = [50, 100, 200, 500]
max_d = [4, 5, 6, 7]
gamma = [0.1, 0.2, 0.3, 0.4]
subco = [0.6, 1.0]
subsa = [0.6, 1.0]
len(list(product(n_est, max_d, gamma, subco, subsa)))
import pandas as pd
df = pd.DataFrame(list(product(n_est, max_d, gamma, subco, subsa)), 
                  columns=['n_estimators','max_depth','gamma','colsample_bytree','subsample'])

df.to_csv('xgboost_config_list.csv',index=None)
