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

#include "RealSlicerOV.h"
//#include "DiscreteSlicer.h"
#include "SliceFactoryOV.h"

#include <sampler/MutableSampler.h>
#include <sampler/MutableSampleMethod.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <vector>

using std::vector;
using std::string;
namespace jags {

bool
SliceFactoryOV::canSample(StochasticNode * node, Graph const &graph) const
{
    return RealSlicerOV::canSample(node);
}

Sampler *SliceFactoryOV::makeSampler(StochasticNode *snode,
				     Graph const &graph) const
{
    unsigned int nchain = snode->nchain();
    vector<MutableSampleMethod*> methods(nchain, 0);

    vector<StochasticNode*> sv;
    sv.push_back(snode);
    GraphView * gv = new GraphView(sv, graph);

    for (unsigned int ch = 0; ch < nchain; ++ch) {

	methods[ch] = new RealSlicerOV(gv, ch);

    }
    string name;
    name = "base::RealSlicer";

    vector<StochasticNode*> nodes(1, snode);
    return new MutableSampler(gv, methods, name);
}

std::string SliceFactoryOV::name() const
{
    return "SliceFactoryOV";
}

}

