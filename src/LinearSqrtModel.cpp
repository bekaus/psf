#include <cmath>
#include <psf/Error.h>
#include "psf/PeakParameter.h"

namespace psf
{

unsigned int LinearSqrtModel::numberOfParameters() {
    return numberOfParameters_;
}

void LinearSqrtModel::setParameter(unsigned index, double value) {
    psf_precondition(index < numberOfParameters(), "LinearSqrtModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double LinearSqrtModel::getParameter(unsigned index) {
    psf_precondition(index < numberOfParameters(), "LinearSqrtModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double LinearSqrtModel::at(const double x) const {
    psf_precondition(x >= 0, "LinearSqrtModel::at(): Parameter x has to be >= 0.");
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
    psf_precondition(index < numberOfParameters(), "LinearSqrtModel::setParameter(): Parameter index out-of-range.");    
    a_ = value;
}
double LinearSqrtOriginModel::getParameter(unsigned index) {
    psf_precondition(index < numberOfParameters(), "LinearSqrtModel::getParameter(): Parameter index out-of-range.");
    return a_;
}

double LinearSqrtOriginModel::at(const double x) const {
    psf_precondition(x >= 0, "LinearSqrtOriginModel::at(): Parameter x has to be >= 0.");
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

} /* namespace psf */
