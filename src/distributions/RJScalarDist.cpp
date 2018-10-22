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
##    Copyright © 2018 Frank A. Schmid,                                                         ##
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

#include "RJScalarDist.h"
#include <RNG.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <cmath>
#include <algorithm>

using std::string;
using std::vector;
using std::log;
using std::min;
using std::max;

namespace jags {

double RJScalarDist::calPlower(double lower,
			      vector<double const*> const &parameters) const
{
    //P(X < lower)
    if (_discrete) {
	return p(lower - 1, parameters, true, false);
    }
    else {
	return p(lower, parameters, true, false);
    }
}

double RJScalarDist::calPupper(double upper,
			      vector<double const*> const &parameters) const
{
    //P(X <= upper)
    return p(upper, parameters, true, false);
}


RJScalarDist::RJScalarDist(string const &name, unsigned int npar,
			 Support support, bool discrete)

    : ScalarDist(name, npar, support),  _support(support), _discrete(discrete),
      _npar(npar)
{
}

double
RJScalarDist::typicalValue(vector<double const *> const &parameters,
			  double const *lower, double const *upper) const
{
    double llimit = l(parameters), ulimit = u(parameters);
    double plower = 0, pupper = 1;

    if (lower) {
	llimit = max(llimit, *lower);
	plower = calPlower(llimit, parameters);
    }

    if (upper) {
	ulimit = min(ulimit, *upper);
	pupper = calPupper(ulimit, parameters);
    }

    double pmed = (plower + pupper)/2;
    double med = q(pmed, parameters, true, false);

    //Calculate the log densities
    double dllimit = d(llimit, parameters, true);
    double dulimit = d(ulimit, parameters, true);
    double dmed = d(med, parameters, true);

    //Pick the median if it has the highest density, otherwise pick
    //a point near to (but not on) the boundary
    if (dmed >= dllimit && dmed >= dulimit) {
	return med;
    }
    else if (dulimit > dllimit) {
	return q(0.1 * plower + 0.9 * pupper, parameters, true, false);
    }
    else {
	return q(0.9 * plower + 0.1 * pupper, parameters, true, false);
    }
}

double
RJScalarDist::logDensity(double x, PDFType type,
                        vector<double const *> const &parameters,
                        double const *lower, double const *upper) const
{
    double loglik =  d(x, parameters, true);

    if (lower || upper) {

	if (lower && x < *lower)
	    return JAGS_NEGINF;
	if (upper && x > *upper)
	    return JAGS_NEGINF;
	if (upper && lower && *upper < *lower)
	    return JAGS_NEGINF;

	//Make adjustment for discrete-valued distributions
	double ll = 0;
	if (lower) {
	    ll = _discrete ? *lower - 1 : *lower;
	}

	/* In theory, we just have to subtract log[P(lower <= X <=
           upper)] from the log likelihood. But we need to work around
           numerical problems. */

	bool have_lower = lower && p(ll, parameters, true, false) > 0;
	bool have_upper = upper && p(*upper, parameters, false, false) > 0;

	if (have_lower && have_upper) {
	    if (p(ll, parameters, false, false) < 0.5) {
		//Use upper tail
		loglik -= log(p(ll, parameters, false, false) -
			      p(*upper, parameters, false, false));
	    }
	    else {
		//Use lower tail
		loglik -= log(p(*upper, parameters, true, false) -
			      p(ll, parameters, true, false));
	    }
	}
	else if (have_lower) {
	    loglik -= p(ll, parameters, false, true);
	}
	else if (have_upper) {
	    loglik -= p(*upper, parameters, true, true);
	}
    }

    return loglik;
}


double
RJScalarDist::randomSample(vector<double const *> const &parameters,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    if (lower || upper) {

	double plower = 0, pupper = 1;
	if (lower) {
	    plower = calPlower(*lower, parameters);
	}
	if (upper) {
	    pupper = calPupper(*upper, parameters);
	}

	double u = plower + rng->uniform() * (pupper - plower);
	return q(u, parameters, true, false);
    }
    else {
	return r(parameters, rng);
    }

}

bool RJScalarDist::canBound() const
{
    return true;
}

bool RJScalarDist::isDiscreteValued(vector<bool> const &mask) const
{
    return _discrete;
}

bool RJScalarDist::discrete() const
{
    return _discrete;
}

unsigned int RJScalarDist::npar() const
{
    return _npar;
}
}
