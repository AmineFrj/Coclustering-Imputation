# -*- coding: utf-8 -*-

"""
The :mod:`coclust.coclustering.coclust_mod` module provides an implementation
of a co-clustering algorithm by direct maximization of graph modularity.
"""

# Author: Francois Role <francois.role@gmail.com>
#         Stanislas Morbieu <stanislas.morbieu@gmail.com>

# License: BSD 3 clause

import numpy as np
from sklearn.utils import check_random_state, check_array

import time
import sys

from ..initialization import random_init
from .base_diagonal_coclust import BaseDiagonalCoclust


class CoclustMod(BaseDiagonalCoclust):
    """Co-clustering by direct maximization of graph modularity.

    Parameters
    ----------
    n_clusters : int, optional, default: 2
        Number of co-clusters to form

    init : numpy array or scipy sparse matrix, \
        shape (n_features, n_clusters), optional, default: None
        Initial column labels

    max_iter : int, optional, default: 20
        Maximum number of iterations

    n_init : int, optional, default: 1
        Number of time the algorithm will be run with different
        initializations. The final results will be the best output of `n_init`
        consecutive runs in terms of modularity.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    tol : float, default: 1e-9
        Relative tolerance with regards to modularity to declare convergence

    Attributes
    ----------
    row_labels_ : array-like, shape (n_rows,)
        Bicluster label of each row

    column_labels_ : array-like, shape (n_cols,)
        Bicluster label of each column

    modularity : float
        Final value of the modularity

    modularities : list
        Record of all computed modularity values for all iterations

    References
    ----------
    * Ailem M., Role F., Nadif M., Co-clustering Document-term Matrices by \
    Direct Maximization of Graph Modularity. CIKM 2015: 1807-1810
    """

    def __init__(self, n_clusters=2, init=None, max_iter=20, n_init=1, missing_indexes=None,
                 tol=1e-9, random_state=None):
        self.n_clusters = n_clusters
        self.init = init
        self.max_iter = max_iter
        self.n_init = n_init
        self.tol = tol
        self.random_state = random_state

        self.row_labels_ = None
        self.column_labels_ = None
        self.modularity = -np.inf
        self.modularities = []

        self.missing = missing_indexes

    def fit(self, X, y=None):
        """Perform co-clustering by direct maximization of graph modularity.

        Parameters
        ----------
        X : numpy array or scipy sparse matrix, shape=(n_samples, n_features)
            Matrix to be analyzed
        """

        random_state = check_random_state(self.random_state)

        check_array(X, accept_sparse=True, dtype="numeric", order=None,
                    copy=False, force_all_finite=True, ensure_2d=True,
                    allow_nd=False, ensure_min_samples=self.n_clusters,
                    ensure_min_features=self.n_clusters,
                    # warn_on_dtype=False, 
                    # FutureWarning: 'warn_on_dtype' is deprecated in version 0.21 
                    # and will be removed in 0.23. Don't set `warn_on_dtype` to remove this warning.
                    estimator=None)

        if type(X) == np.ndarray:
            X = np.matrix(X)

        X = X.astype(float)

        modularity = self.modularity
        modularities = []
        row_labels_ = None
        column_labels_ = None
        self.X = X

        seeds = random_state.randint(np.iinfo(np.int32).max, size=self.n_init)
        for seed in seeds:
            self._fit_single(self.X, seed, y)
            if np.isnan(self.modularity):
                raise ValueError("matrix may contain unexpected NaN values")
            # remember attributes corresponding to the best modularity
            if (self.modularity > modularity):
                modularity = self.modularity
                modularities = self.modularities
                row_labels_ = self.row_labels_
                column_labels_ = self.column_labels_

        # update attributes
        self.modularity = modularity
        self.modularities = modularities
        self.row_labels_ = row_labels_
        self.column_labels_ = column_labels_

        if (self.missing is not None):
            sys.stdout.write('\rProgress: ['+'='*16+']   ## With Imputation ##\n')
        else:
            sys.stdout.write('\rProgress: ['+'='*16+']   ## Without Imputation ##\n')

        return self

    def _fit_single(self, X, random_state, y=None):
        """Perform one run of co-clustering by direct maximization of graph
        modularity.

        Parameters
        ----------
        X : numpy array or scipy sparse matrix, shape=(n_samples, n_features)
            Matrix to be analyzed
        """

        if self.init is None:
            W = random_init(self.n_clusters, X.shape[1], random_state)
        else:
            W = np.matrix(self.init, dtype=float)

        Z = np.zeros((X.shape[0], self.n_clusters))

        # Compute the modularity matrix
        row_sums = np.matrix(self.X.sum(axis=1))
        col_sums = np.matrix(self.X.sum(axis=0))
        N = float(self.X.sum())
        indep = (row_sums.dot(col_sums)) / N

        # B is a numpy matrix
        B = self.X - indep

        self.modularities = []

        # Loop
        m_begin = float("-inf")
        change = True
        iteration = 0
        while change:
            change = False

            # Reassign rows
            BW = B.dot(W)
            for idx, k in enumerate(np.argmax(BW, axis=1)):
                Z[idx, :] = 0
                Z[idx, k] = 1
            # Imputation
            prog = int(iteration*16/self.max_iter)
            if (self.missing is not None):
                sys.stdout.write('\rProgress: ['+'='*prog+'>'+' '*(15-prog)+']   ## With Imputation ##')
                self._impute(Z=Z, W=W)
            else:
                sys.stdout.write('\rProgress: ['+'='*prog+'>'+' '*(15-prog)+']   ## Without Imputation ##')

            # Reassign columns
            BtZ = (B.T).dot(Z)
            for idx, k in enumerate(np.argmax(BtZ, axis=1)):
                W[idx, :] = 0
                W[idx, k] = 1
            # Imputation
            if (self.missing is not None):
                self._impute(Z=Z, W=W)

            k_times_k = (Z.T).dot(BW)
            m_end = np.trace(k_times_k)
            iteration += 1
            if (np.abs(m_end - m_begin) > self.tol and
                    iteration < self.max_iter):
                self.modularities.append(m_end/N)
                m_begin = m_end
                change = True

        self.row_labels_ = np.argmax(Z, axis=1).tolist()
        self.column_labels_ = np.argmax(W, axis=1).tolist()
        self.btz = BtZ
        self.bw = BW
        self.modularity = m_end / N
        self.nb_iterations = iteration

    def _impute(self, Z, W):
        """Imputes missing values from given indexes

        Parameters
        ----------
        
        """
        A = np.zeros((self.n_clusters,self.n_clusters))
        for i in range(self.n_clusters):
            Z_ind = np.where(Z[:,i] == 1)
            for j in range(self.n_clusters):
                W_ind = np.where(W[:,j] == 1)
                row_idx = Z_ind[0]
                col_idx = W_ind[0]
                tmp = self.X[row_idx[:, None], col_idx]
                if len(tmp) and tmp.sum() > 0:
                    A[i,j] = tmp.mean()
                else :
                    A[i,j] = 0
        A = np.nan_to_num(A)
        # To modify X values use lil_matrix much efficient
        for i in self.missing.values:
            cluster_row = Z[i[0],:].tolist().index(1)
            cluster_col = W[i[1],:].tolist().index(1)
            self.X[i[0],i[1]] = A[cluster_row, cluster_col]

    def get_assignment_matrix(self, kind, i):
        """Returns the indices of 'best' i cols of an assignment matrix
        (row or column).

        Parameters
        ----------
        kind : string
             Assignment matrix to be used: rows or cols

        Returns
        -------
        numpy array or scipy sparse matrix
            Matrix containing the i 'best' columns of a row or column
            assignment matrix
        """
        if kind == "rows":
            s_bw = np.argsort(self.bw)
            return s_bw[:, -1:-(i+1):-1]
        if kind == "cols":
            s_btz = np.argsort(self.btz)
            return s_btz[:, -1:-(i+1):-1]
