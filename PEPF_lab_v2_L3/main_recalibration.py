"""
Main script to run the recalibration experiments
"""
# Author: Alessandro Brusaferri
# License: Apache-2.0 license

import os
import pandas as pd
import numpy as np
os.environ["TF_USE_LEGACY_KERAS"]="1"
from tools.PrTSF_Recalib_tools import PrTsfRecalibEngine, load_data_model_configs
from tools.prediction_quantiles_tools import plot_quantiles

#--------------------------------------------------------------------------------------------------------------------
def compute_pinball_scores(y_true, pred_quantiles, quantiles_levels):
    """
    Utility function to compute the pinball score on the test results
    return: pinball scores computed for each quantile level and each step in the pred horizon
    """
    score = []
    for i, q in enumerate(quantiles_levels):
        error = np.subtract(y_true, pred_quantiles[:, :, i])
        loss_q = np.maximum(q * error, (q - 1) * error)
        score.append(np.expand_dims(loss_q,-1))
    score = np.mean(np.concatenate(score, axis=-1), axis=0)
    return score

#--------------------------------------------------------------------------------------------------------------------

# compute the winkler score
def compute_winkler_scores(y_true, pred_quantiles, quantiles_levels):
    """
    Utility function to compute the winkler score on the test results
    return: winkler scores computed for each quantile level and each step in the pred horizon
    """
    score = []
    # loop over only half the quantiles
    # the other half is symmetric, winkler is only applied to quantiles above (or below) the median
    for i, tau in enumerate(quantiles_levels[:len(quantiles_levels)//2]):
        # get the upper and lower quantiles
        L_tau = pred_quantiles[:, :, i]
        U_tau = pred_quantiles[:, :, -i-1]
        # compute the quantile width
        delta_tau = np.subtract(U_tau, L_tau)
        # compute the errors
        error_L = np.subtract(L_tau, y_true)
        error_U = np.subtract(y_true, U_tau)
        # compute the winkler score
        # use tau, not 1-tau, as the quantiles are symmetric
        loss_q = delta_tau + 2 / tau * (
            np.maximum(error_L, np.zeros(error_L.shape))
            + np.maximum(error_U, np.zeros(error_U.shape))
        )
        score.append(np.expand_dims(loss_q,-1))
    score = np.mean(np.concatenate(score, axis=-1), axis=0)

    return score

#--------------------------------------------------------------------------------------------------------------------
# Set PEPF task to execute
PF_task_name = 'EM_price'
# Set Model setup to execute
exper_setup = 'QR-DNN-ARCSINH'

#---------------------------------------------------------------------------------------------------------------------
# Set run configs

# fix the seed for reproducibility
np.random.seed(42)

# run_id = 'recalib_opt_grid_1_1'
run_id = 'recalib_opt_random_1_2'
# Load hyperparams from file (select: load_tuned or optuna_tuner)
hyper_mode = 'optuna_tuner'
# Plot train history flag
plot_train_history=False
plot_weights=False
# Apply sinh transformation to the target variable
apply_arcsinh_transf = False

#---------------------------------------------------------------------------------------------------------------------
# Load experiments configuration from json file
configs=load_data_model_configs(task_name=PF_task_name, exper_setup=exper_setup, run_id=run_id)

# Load dataset
# ***: This lets us change the dataset before passing it directly to the recalibration engine
# ***: e.g. preprocessing, feature selection, etc.
dir_path = os.getcwd()
ds = pd.read_csv(os.path.join(dir_path, 'data', 'datasets', configs['data_config'].dataset_name))
ds.set_index(ds.columns[0], inplace=True)

if apply_arcsinh_transf:
    ds['TARG__'+PF_task_name] = np.arcsinh(ds['TARG__'+PF_task_name])

#---------------------------------------------------------------------------------------------------------------------
# Instantiate recalibratione engine
PrTSF_eng = PrTsfRecalibEngine(dataset=ds,
                               data_configs=configs['data_config'],
                               model_configs=configs['model_config'])

# Get model hyperparameters (previously saved or by tuning)
model_hyperparams = PrTSF_eng.get_model_hyperparams(method=hyper_mode, optuna_m=configs['model_config']['optuna_m'])

# Exec recalib loop over the test_set samples, using the tuned hyperparams
test_predictions = PrTSF_eng.run_recalibration(model_hyperparams=model_hyperparams,
                                               plot_history=plot_train_history,
                                               plot_weights=plot_weights)

# apply inverse sinh transformation to all the predictions
if apply_arcsinh_transf:
    test_predictions = np.sinh(test_predictions)

#--------------------------------------------------------------------------------------------------------------------
# Compute pinball score
quantiles_levels = PrTSF_eng.model_configs['target_quantiles']
pred_steps = configs['model_config']['pred_horiz']

pinball_scores = compute_pinball_scores(y_true=test_predictions[PF_task_name].to_numpy().reshape(-1,pred_steps),
                                        pred_quantiles=test_predictions.loc[:,test_predictions.columns != PF_task_name].
                                        to_numpy().reshape(-1, pred_steps, len(quantiles_levels)),
                                        quantiles_levels=quantiles_levels)

# print the Pinball score as a table
pinbal_df = pd.DataFrame(pinball_scores, columns=[f'q_{q}' for q in quantiles_levels], index=[f'Hour {i+1}' for i in range(pred_steps)])
print('--- Pinball Scores ---')
print(pinbal_df)

# save pinball scores to file
pinbal_df.to_csv(os.path.join(dir_path, 'experiments', 'tasks', PF_task_name, exper_setup, run_id, 'pinball_scores.csv'))

# Compute winkler score
winkler_scores = compute_winkler_scores(y_true=test_predictions[PF_task_name].to_numpy().reshape(-1,pred_steps),
                                        pred_quantiles=test_predictions.loc[:,test_predictions.columns != PF_task_name].
                                        to_numpy().reshape(-1, pred_steps, len(quantiles_levels)),
                                        quantiles_levels=quantiles_levels)

# print the Winkler score as a table
winkler_df = pd.DataFrame(winkler_scores, columns=[f'q_{q}' for q in quantiles_levels[:len(quantiles_levels)//2]],
    index=[f'Hour {i+1}' for i in range(pred_steps)])
print('--- Winkler Scores ---')
print(winkler_df)

# save winkler scores to file
winkler_df.to_csv(os.path.join(dir_path, 'experiments', 'tasks', PF_task_name, exper_setup, run_id, 'winkler_scores.csv'))

#--------------------------------------------------------------------------------------------------------------------
# Plot test predictions
plot_quantiles(test_predictions, target=PF_task_name)

# stop the execution until enter is pressed
input("Press Enter to continue...")

#--------------------------------------------------------------------------------------------------------------------
print('Done!')
