"""
Ensemble model
"""

# Author: Alessandro Brusaferri
# License: Apache-2.0 license

import sys
import numpy as np
import tensorflow as tf
#import tensorflow_probability as tfp
#from tensorflow_probability import distributions as tfd
import matplotlib.pyplot as plt

from tools.models.ARX import ARXRegressor

def get_model_class_from_conf(conf):
    """
    Map the model class depending on the config name
    """
    if conf == 'ARX':
        model_class = ARXRegressor
    else:
        sys.exit('ERROR: unknown model_class')
    return model_class


def regression_model(settings, sample_x):
    """
    Wrapper to the regression model
    :param settings: model configurations, sample_x: input sample to derive the model input shape (first dimension has to be 1)
    :param sample_x: input sample, to derive the model input shape (first dimension has to be 1)
    :return: instantiated model
    """
    # Currently direct link to TF, Future dev pytorch
    return TensorflowRegressor(settings=settings, sample_x=sample_x)


class TensorflowRegressor():
    """
    Implementation of the Tenforflow regressor
    """
    def __init__(self, settings, sample_x):
        self.settings = settings
        self.x_columns_names = settings['x_columns_names']
        self.pred_horiz = settings['pred_horiz']

        tf.keras.backend.clear_session()
        # Map the loss to be used
        if settings['PF_method']=='point':
            loss = 'mae'
        else:
            sys.exit('ERROR: unknown PF_method config!')

        # Instantiate the model
        if  settings['model_class']=='ARX':
            # get input size for the chosen model architecture
            settings['input_size']=ARXRegressor.build_model_input_from_series(x=sample_x,
                col_names=self.x_columns_names, pred_horiz=self.pred_horiz).shape[1]
            # Build the model architecture
            self.regressor = ARXRegressor(settings, loss)
        else:
            sys.exit('ERROR: unknown model_class')

        # Map handler to convert distributional output to quantiles or distribution parameters
        self.output_handler =self.__quantiles_out__

    def fit(self, train_x, train_y, val_x, val_y, verbose=0, pruning_call=None, plot_history=False):
        history = self.regressor.fit(train_x, train_y, val_x, val_y, verbose=0, pruning_call=None)
        if plot_history:
            plt.plot(history.history['loss'], label='train_loss')
            plt.plot(history.history['val_loss'], label='vali_loss')
            plt.grid()
            plt.legend()
            # save the figure with a random name
            # plt.savefig('train_history' + str(np.random.randint(0, 10000)) + '.png')
            plt.show()

    def predict(self, x):
        return self.output_handler(self.regressor.predict(x))

    def evaluate(self, x, y):
        return self.regressor.evaluate(x=x, y=y)

    def __quantiles_out__(self, preds):
        # Expand dimension to enable concat in ensemble
        return tf.expand_dims(preds, axis=2)

    def plot_weights(self):
        self.regressor.plot_weights()

class Ensemble():
    """
    Tensorflow ensemble wrapper
    """
    def __init__(self, settings):
        # store configs for internal use
        self.settings = settings
        # map the methods to use for aggretation and quantile building depending on the configs
        if (self.settings['PF_method'] == 'point'):
            self.ensemble_aggregator = self.__aggregate_de_quantiles__
            self._build_test_PIs = self.__get_qr_PIs__
        else:
            sys.exit('ERROR: Ensemble config not supported!')

    def aggregate_preds(self, ens_comp_preds):
        # link function to the specific aggregator
        return self.ensemble_aggregator(ens_comp_preds=ens_comp_preds)

    def get_preds_test_quantiles(self, preds_test):
        # link function to the specific PI builder
        return self._build_test_PIs(preds_test=preds_test, settings=self.settings)

    @staticmethod
    def __aggregate_de__(ens_comp_preds):
        # aggregate by concatenation, for point a distributional settings
        return np.concatenate(ens_comp_preds, axis=2)

    @staticmethod
    def __aggregate_de_quantiles__(ens_comp_preds):
        # aggregate by a uniform vincentization
        return np.mean(np.concatenate(ens_comp_preds, axis=2), axis=2)

    @staticmethod
    def __get_qr_PIs__(preds_test, settings):
        # simply flatten in temporal dimension
        return preds_test.reshape(-1, preds_test.shape[-1])