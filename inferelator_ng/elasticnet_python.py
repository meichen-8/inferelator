import numpy as np

from inferelator_ng import utils
from inferelator_ng import regression
from inferelator_ng import tfa_workflow
from inferelator_ng.kvs_controller import ownCheck
from sklearn.linear_model import ElasticNetCV

ELASTICNET_PARAMETERS = dict(l1_ratio=[0.5, 0.7, 0.9],
                             eps=0.001,
                             n_alphas=50,
                             alphas=None,
                             fit_intercept=True,
                             normalize=False,
                             precompute='auto',
                             max_iter=1000,
                             tol=0.001,
                             cv=3,
                             copy_X=True,
                             verbose=0,
                             n_jobs=1,
                             positive=False,
                             random_state=99,
                             selection='random')

MIN_COEF = 0.1


def elastic_net(X, Y, params):
    """

    :param X: np.ndarray [K x N]
    :param Y: np.ndarray [1 x N]
    :param params: dict
    :return:
    """
    (K, N) = X.shape
    X = X.T  # Make X into [N, K]
    Y = Y.flatten()  # Make Y into [N, ]

    # Fit the linear model using the elastic net
    model = ElasticNetCV(**params).fit(X, Y)

    # Set coefficients below threshold to 0
    coefs = model.coef_  # Get all model coefficients [K, ]
    coefs[np.abs(coefs) < MIN_COEF] = 0.  # Threshold coefficients
    coef_nonzero = coefs != 0  # Create a boolean array where coefficients are nonzero [K, ]

    # If there are non-zero coefficients, redo the linear regression with them alone
    # And calculate beta_resc
    if coef_nonzero.sum() > 0:
        utils.make_array_2d(Y)
        return regression.recalculate_betas_from_selected(X, Y, coef_nonzero)
    else:
        return dict(pp=np.repeat(True, K).tolist(),
                    betas=np.zeros(K),
                    betas_resc=np.zeros(K))


class ElasticNet(regression.BaseRegression):
    params = ELASTICNET_PARAMETERS

    def __init__(self, X, Y, kvs, chunk=25):
        """
        Create a ElasticNet object for regularized regression

        :param X: pd.DataFrame [K x N]
            Expression / Activity data
        :param Y: pd.DataFrame [G x N]
            Response data
        """
        self.kvs = kvs
        self.chunk = chunk

        # Get the IDs and total count for the genes and predictors
        self.K = X.shape[0]
        self.tfs = X.index.values.tolist()
        self.G = Y.shape[0]
        self.genes = Y.index.values.tolist()

        # Rescale input data
        self.X = self._scale(X)
        utils.Debug.vprint("Predictor matrix {} ready".format(X.shape))
        self.Y = self._scale(Y)
        utils.Debug.vprint("Response matrix {} ready".format(Y.shape))

    def run(self):
        """
        Execute regression separately on each response variable in the data

        :return: pd.DataFrame [G x K], pd.DataFrame [G x K]
            Returns the regression betas and beta error reductions for all threads if this is the master thread (rank 0)
            Returns None, None if it's a subordinate thread
        """
        regression_data = []

        # For every response variable G, check to see if this thread should run BBSR for that variable
        # If it should (ownCheck is TRUE), slice the data for that response variable
        # And send the values (as an ndarray because pd.df indexing is glacially slow) to bayes_stats.bbsr
        # Keep a list of the resulting regression results
        oc = ownCheck(self.kvs, self.kvs.rank, chunk=25)
        for j in range(self.G):
            if next(oc):
                level = 0 if j % 100 == 0 else 2
                utils.Debug.vprint("Regression on {gn} [{i} / {total}]".format(gn=self.Y.index[j],
                                                                               i=j,
                                                                               total=self.G), level=level)
                data = elastic_net(self.X.values, self.Y.ix[j, :].values, self.params)
                data['ind'] = j
                regression_data.append(data)

        # Put the regression results that this thread has calculated into KVS
        self.kvs.put('plist', (self.kvs.rank, regression_data))
        self.kvs.sync_processes("bootstrap")

        # If this is the master thread, pile the regression betas into dataframes and return them
        if self.kvs.is_master:
            return self.pileup_data()
        else:
            return None, None


class ElasticNetRunner:
    def run(self, X, Y, kvs):
        return ElasticNet(X, Y, kvs).run()


class MEN_Workflow(tfa_workflow.TFAWorkFlow):
    # Drivers
    regression_driver = ElasticNetRunner

    def run_bootstrap(self, bootstrap):
        X = self.design.iloc[:, bootstrap]
        Y = self.response.iloc[:, bootstrap]
        utils.Debug.vprint('Calculating betas using MEN', level=0)
        self.kvs.sync_processes("pre-bootstrap")
        return self.regression_driver().run(X, Y, self.kvs)
