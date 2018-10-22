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


//#include <config.h>
#include "DTOV.h"

#include <cmath>

#include <JRmath.h>

using std::vector;

#define MU(par) (*par[0])
#define TAU(par) (*par[1])
#define DF(par) (*par[2])

namespace jags {

DTOV::DTOV()
    : RJScalarDist("dtOV", 3, DIST_UNBOUNDED)
{}

bool DTOV::checkParameterValue (vector<double const *> const &par) const
{
    return (TAU(par) > 0 && DF(par) > 0);
}

double DTOV::d(double x, vector<double const *> const &par, bool give_log) const
{
    x = (x - MU(par)) * sqrt(TAU(par));
    if (give_log) {
	return dt(x, DF(par), 1) + log(TAU(par))/2;
    }
    else {
	return dt(x, DF(par), 0) * sqrt(TAU(par));
    }
}

double DTOV::p(double x, vector<double const *> const &par, bool lower,
	     bool use_log) const
{
    return pt((x - MU(par)) * sqrt(TAU(par)), DF(par), lower, use_log);
}

double DTOV::q(double p, vector<double const *> const &par, bool lower,
	     bool log_p) const
{
    return MU(par) + qt(p, DF(par), lower, log_p) / sqrt(TAU(par));
}

double DTOV::r(vector<double const *> const &par, RNG *rng) const
{
    return rt(DF(par), rng) / sqrt(TAU(par)) + MU(par);
}
}
