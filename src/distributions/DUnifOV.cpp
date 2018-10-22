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
#include "DUnifOV.h"

#include <cmath>
#include <RNG.h>

using std::vector;
using std::log;

#define LOWER(par) (*par[0])
#define UPPER(par) (*par[1])

namespace jags {

DUnifOV::DUnifOV()
    : ScalarDist("dunifOV", 2, DIST_SPECIAL)
{}

bool  DUnifOV::checkParameterValue (vector<double const *> const &par) const
{
    return (LOWER(par) < UPPER(par));
}

double DUnifOV::logDensity(double x, PDFType type,
                           vector<double const *> const &par,
                           double const *lower, double const *upper) const
{
    return log(UPPER(par) - LOWER(par));
}

double DUnifOV::randomSample(vector<double const *> const &par,
			   double const *lower, double const *upper,
			   RNG *rng) const
{
    return LOWER(par) + rng->uniform() * (UPPER(par) - LOWER(par));
}

double DUnifOV::typicalValue(vector<double const *> const &par,
			   double const *lower, double const *upper) const
{
    return (LOWER(par) + UPPER(par))/2;
}

double DUnifOV::l(vector<double const*> const &par) const
{
    return LOWER(par);
}

double DUnifOV::u(vector<double const*> const &par) const
{
    return UPPER(par);
}

bool DUnifOV::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[0] && fixmask[1]; //Lower and upper bounds fixed
}
}
