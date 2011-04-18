/* $Id: SqrtModel.cpp 2253 2009-08-19 11:43:25Z bkausler $ */

/*
 * SqrtModel.cpp
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

#include <ms++/PeakParameter.h>

using namespace ms;

unsigned int SqrtModel::numberOfParameters() {
    return numberOfParameters_;
}

void SqrtModel::setParameter(unsigned index, double value) {
    mspp_precondition(index < numberOfParameters(), "SqrtModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double SqrtModel::getParameter(unsigned index) {
    mspp_precondition(index < numberOfParameters(), "SqrtModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double SqrtModel::at(const double x) const {
    mspp_precondition(x >= 0, "SqrtModel::at(): Parameter x hast to be >= 0.");
    return a_ * std::sqrt(x) + b_;
}

GeneralizedSlope SqrtModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {std::sqrt(x), 1., 0.};
    return GeneralizedSlope(slope, slope + 3);
}

// setter / getter
void SqrtModel::setA(const double a) {
    a_ = a;
}
double SqrtModel::getA() const {
    return a_;
}

void SqrtModel::setB(const double b) {
    b_ = b;
}
double SqrtModel::getB() const {
    return b_;
}
