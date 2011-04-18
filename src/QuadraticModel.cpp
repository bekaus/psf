/* $Id: QuadraticModel.cpp 2253 2009-08-19 11:43:25Z bkausler $ */

/*
 * QuadraticModel.cpp
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

#include <cmath>

#include <ms++/PeakParameter.h>

using namespace ms;

unsigned int QuadraticModel::numberOfParameters() {
    return numberOfParameters_;
}

void QuadraticModel::setParameter(unsigned index, double value) {
    mspp_precondition(index < numberOfParameters(), "QuadraticModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double QuadraticModel::getParameter(unsigned index) {
    mspp_precondition(index < numberOfParameters(), "QuadraticModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double QuadraticModel::at(const double x) const {
    return a_ * x*x + b_;
}

GeneralizedSlope QuadraticModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {x * x, 1., 0.};
    return GeneralizedSlope(slope, slope + 3);
}

// setter / getter
void QuadraticModel::setA(const double a) {
    a_ = a;
}
double QuadraticModel::getA() const {
    return a_;
}

void QuadraticModel::setB(const double b) {
    b_ = b;
}
double QuadraticModel::getB() const {
    return b_;
}
