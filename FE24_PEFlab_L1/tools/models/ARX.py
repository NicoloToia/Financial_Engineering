"""
ARX model class
"""
import matplotlib.pyplot as plt

# Author: Alessandro Brusaferri
# License: Apache-2.0 license

from tools.data_utils import features_keys
import numpy as np
import tensorflow as tf
from typing import List
import os
import shutil


class ARXRegressor:
    def __init__(self, settings, loss):
        self.settings = settings
        self.__build_model__(loss)

    def __build_model__(self, loss):
        print('to be implemented')

    def fit(self, train_x, train_y, val_x, val_y, verbose=0, pruning_call=None):
        print('to be implemented')

    def predict(self, x):
        print('to be implemented')

    def evaluate(self, x, y):
        print('to be implemented')

    @staticmethod
    def build_model_input_from_series(x, col_names: List, pred_horiz: int):
        print('to be implemented')

    @staticmethod
    def get_hyperparams_trial(trial, settings):
        settings['l1'] = trial.suggest_float('l1', 1e-7, 1e-1)
        settings['lr'] = trial.suggest_float('lr', 1e-5, 1e-1, log=True)
        return settings

    @staticmethod
    def get_hyperparams_searchspace():
        return {'l1': [1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
                'lr': [1e-4, 1e-3, 1e-2]}

    @staticmethod
    def get_hyperparams_dict_from_configs(configs):
        model_hyperparams = {
            'l1': configs['l1'],
            'lr': configs['lr']
        }
        return model_hyperparams

    def plot_weights(self):
        print('To be implemented!')


