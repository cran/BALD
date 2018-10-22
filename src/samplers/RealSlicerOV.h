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


#ifndef REAL_SLICER_OV_H_
#define REAL_SLICER_OV_H_

#include <sampler/MutableSampleMethod.h>
#include <sampler/GraphView.h>
namespace jags {
class StochasticNode;


/**
 * Slice sampler for real-valued distributions
 */
class RealSlicerOV : public MutableSampleMethod
{
    GraphView const *_gv;
    unsigned int _chain;

    double _probOfOverrelaxed;
    unsigned int _endpointAccuracy;

    unsigned int _maxOV;
    double _widthOV;
    unsigned int _iterOV;
    double _sumdiffOV;

    bool _adaptOV;



public:
    /**
     * Constructor for Slice Sampler
     * @param node Node to sample
     * @param width Initial width of slice
     * @param maxwidth Maximal width of slice as a multiple of the width
     * parameter
     * @param nburn Length of burnin
     */
    RealSlicerOV(GraphView const * gv, unsigned int chain);
    double value() const;
    void setValue(double value) const;
    void getLimits(double *lower, double *upper) const;
    virtual void update(RNG *rng);
    void updateOverrelaxed(RNG *rng);
    void updateStep(RNG *rng);
    virtual std::string name() const;
    static bool canSample(StochasticNode const *node);
    virtual bool isAdaptive() const;
    virtual void adaptOff();
    virtual bool checkAdaptation() const;
    double logFullConditional() const;
};

}

#endif /* REAL_SLICER_OV_H_ */

