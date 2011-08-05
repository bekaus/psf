/* $Id: ConstantModel.cpp 2253 2009-08-19 11:43:25Z bkausler $ */

/*
 * ConstantModel.cpp
 *
 * Copyright (c) 2009 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 *
 * This file is part of ms++.
 *
 * ms++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ms++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ms++. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <ms++/Error.h>
#include "psf/PeakParameter.h"

namespace psf
{
unsigned int ConstantModel::numberOfParameters() {
    return numberOfParameters_;
}

void ConstantModel::setParameter(unsigned index, double value) {
    mspp_precondition(index < numberOfParameters(), "ConstantModel::setParameter(): Parameter index out-of-range.");
    a_ = value;
}
double ConstantModel::getParameter(unsigned index) {
    mspp_precondition(index < numberOfParameters(), "ConstantModel::getParameter(): Parameter index out-of-range.");
    return a_;
}

double ConstantModel::at(const double x) const {
    return a_;
}

GeneralizedSlope ConstantModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {1., 0.};
    return GeneralizedSlope(slope, slope + 2);
}

// setter / getter
void ConstantModel::setA(const double a) {
    a_ = a;
}
double ConstantModel::getA() const {
    return a_;
}

} /* namespace psf */
