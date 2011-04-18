/* $Id: LinearSqrtModel.cpp 2253 2009-08-19 11:43:25Z bkausler $ */

/*
 * LinearSqrtModel.cpp
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

#include <ms++/Error.h>

#include "psf/PeakParameter.h"

namespace ms
{

unsigned int LinearSqrtModel::numberOfParameters() {
    return numberOfParameters_;
}

void LinearSqrtModel::setParameter(unsigned index, double value) {
    mspp_precondition(index < numberOfParameters(), "LinearSqrtModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double LinearSqrtModel::getParameter(unsigned index) {
    mspp_precondition(index < numberOfParameters(), "LinearSqrtModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double LinearSqrtModel::at(const double x) const {
    mspp_precondition(x >= 0, "LinearSqrtModel::at(): Parameter x has to be >= 0.");
    return a_ * x * std::sqrt(x) + b_;
}

GeneralizedSlope LinearSqrtModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {x * std::sqrt(x), 1., 0.};
    return GeneralizedSlope(slope, slope + 3);
}

// setter / getter
void LinearSqrtModel::setA(const double a) {
    a_ = a;
}
double LinearSqrtModel::getA() const {
    return a_;
}

void LinearSqrtModel::setB(const double b) {
    b_ = b;
}
double LinearSqrtModel::getB() const {
    return b_;
}




unsigned int LinearSqrtOriginModel::numberOfParameters() {
    return numberOfParameters_;
}

void LinearSqrtOriginModel::setParameter(unsigned index, double value) {
    mspp_precondition(index < numberOfParameters(), "LinearSqrtModel::setParameter(): Parameter index out-of-range.");    
    a_ = value;
}
double LinearSqrtOriginModel::getParameter(unsigned index) {
    mspp_precondition(index < numberOfParameters(), "LinearSqrtModel::getParameter(): Parameter index out-of-range.");
    return a_;
}

double LinearSqrtOriginModel::at(const double x) const {
    mspp_precondition(x >= 0, "LinearSqrtOriginModel::at(): Parameter x has to be >= 0.");
    return a_ * x * std::sqrt(x);
}

GeneralizedSlope LinearSqrtOriginModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {x * std::sqrt(x), 0.};
    return GeneralizedSlope(slope, slope + 2);
}

// setter / getter
void LinearSqrtOriginModel::setA(const double a) {
    a_ = a;
}
double LinearSqrtOriginModel::getA() const {
    return a_;
}

} /* namespace ms */
