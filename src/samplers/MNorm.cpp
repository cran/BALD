/*
##################################################################################################
##                                                                                              ##
##    BALD is an R-package.                                                                     ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of development years.     ##
##    It is an update for the older package lossDev as it has been stopped supported.           ##
##                                                                                              ##
##    Copyright � 2018 Frank A. Schmid,                                                         ##
##                                                                                              ##
##    This file is part of BALD.                                                                ##
##                                                                                              ##
##    lossDev is free software: you can redistribute it and/or modify                           ##
##    it under the terms of the GNU General Public License as published by                      ##
##    the Free Software Foundation, either version 3 of the License, or                         ##
##    (at your option) any later version.                                                       ##
##                                                                                              ##
##    This program is distributed in the hope that it will be useful,                           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of                            ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             ##
##    GNU General Public License for more details.                                              ##
##                                                                                              ##
##    You should have received a copy of the GNU General Public License                         ##
##    along with this program.  If not, see <https://www.gnu.org/licenses/>.                    ##
##                                                                                              ##
##################################################################################################
*/

#include <R_ext/Lapack.h>
#include <float.h> //DBL_EPSILON
#include <JRmath.h>
#include <module/ModuleError.h>

#include <vector>
#include <string>

namespace jags {
class RNG;

using std::vector;


static double logdet(double const *a, int n)
{
    // Log determinant of n x n symmetric positive matrix a */

    int N = n*n;
    double *acopy = new double[N];
    for (int i = 0; i < N; i++) {
	acopy[i] = a[i];
    }

    double *w = new double[n];
    int lwork = -1;
    double worktest = 0;
    int info = 0;
    F77_NAME(dsyev)("N","U", &n, acopy, &n, w, &worktest, &lwork, &info);
    if (info != 0) {
	delete [] acopy;
	delete [] w;
	throwRuntimeError("unable to calculate workspace size for dsyev");
    }
    lwork = static_cast<int>(worktest);
    double *work = new double[lwork];
    F77_NAME(dsyev)("N","U", &n, acopy, &n, w, work, &lwork, &info);
    delete [] acopy;
    delete [] work;
    if (info != 0) {
	delete [] w;
	throwRuntimeError("unable to calculate eigenvalues in dsyev");
    }

    if (w[0] <= 0) {
	throwRuntimeError("Non positive definite matrix in call to logdet");
    }

    double logdet = 0;
    for (int i = 0; i < n; i++) {
	logdet += log(w[i]);
    }
    delete [] w;

    return logdet;
}

double MNorm_logLikelihood(double const *x, unsigned int m,
			    vector<double const *> const &parameters,
			    vector<vector<unsigned int> > const &dims,
			    double const *lower, double const *upper)
{
    double const * mu = parameters[0];
    double const * T = parameters[1];

    double loglik = logdet(T, m)/2;
    double * delta = new double[m];
    for (unsigned int i = 0; i < m; ++i) {
	delta[i] = x[i] - mu[i];
	loglik -= (delta[i] * T[i + i * m] * delta[i])/2;
	for (unsigned int j = 0; j < i; ++j) {
	    loglik -= (delta[i] * T[i + j * m] * delta[j]);
	}
    }
    delete [] delta;

    return loglik;
}

void MNorm_randomsample(double *x, double const *mu, double const *T,
			  bool prec, int nrow, RNG *rng)
{
    int N = nrow*nrow;
    double * Tcopy = new double[N];
    for (int i = 0; i < N; ++i) {
	Tcopy[i] = T[i];
    }
    double * w = new double[nrow];

    int info = 0;
    double worktest;
    int lwork = -1;
    // Workspace query
    F77_NAME(dsyev) ("V", "L", &nrow, Tcopy, &nrow, w, &worktest, &lwork,
		     &info);
    // Now get eigenvalues/vectors with optimal work space
    lwork = static_cast<int>(worktest + DBL_EPSILON);
    double * work = new double[lwork];
    F77_NAME(dsyev) ("V", "L", &nrow, Tcopy, &nrow, w, work, &lwork, &info);
    delete [] work;

    /* Generate independent random normal variates, scaled by
       the eigen values. We reuse the array w. */
    if (prec) {
	for (int i = 0; i < nrow; ++i) {
	    w[i] = rnorm(0, 1/sqrt(w[i]), rng);
	}
    }
    else {
	for (int i = 0; i < nrow; ++i) {
	    w[i] = rnorm(0, sqrt(w[i]), rng);
	}
    }

    /* Now transform them to dependant variates
       (On exit from DSYEV, Tcopy contains the eigenvectors)
    */
    for (int i = 0; i < nrow; ++i) {
	x[i] = mu ? mu[i] : 0;
	for (int j = 0; j < nrow; ++j) {
	    x[i] += Tcopy[i + j * nrow] * w[j];
	}
    }
    delete [] w;
    delete [] Tcopy;
}
}
