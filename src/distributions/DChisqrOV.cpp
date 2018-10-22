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

//#include <config.h>
#include "DChisqrOV.h"

#include <JRmath.h>

using std::vector;

#define DF(par) (*par[0])

namespace jags {

DChisqrOV::DChisqrOV()
    : RJScalarDist("dchisqrOV", 1, DIST_POSITIVE)
{}


bool
DChisqrOV::checkParameterValue (vector<double const *> const &par) const
{
    return (DF(par) > 0);
}

double
DChisqrOV::d(double x, vector<double const *> const &par, bool give_log) const
{
    return dchisq(x, DF(par), give_log);
}

double
DChisqrOV::p(double q, vector<double const *> const &par, bool lower, bool log_p)
  const
{
    return pchisq(q, DF(par), lower, log_p);
}

double
DChisqrOV::q(double p, vector<double const *> const &par, bool lower, bool log_p)
const
{
    return qchisq(p, DF(par), lower, log_p);
}

double DChisqrOV::r(vector<double const *> const &par, RNG *rng) const
{
    return rchisq(DF(par), rng);
}
}
