#include <cmath>

#include <ms++/Error.h>

#include "psf/PeakParameter.h"

using namespace psf;

unsigned int SqrtModel::numberOfParameters() {
    return numberOfParameters_;
}

void SqrtModel::setParameter(unsigned index, double value) {
    psf_precondition(index < numberOfParameters(), "SqrtModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double SqrtModel::getParameter(unsigned index) {
    psf_precondition(index < numberOfParameters(), "SqrtModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double SqrtModel::at(const double x) const {
    psf_precondition(x >= 0, "SqrtModel::at(): Parameter x hast to be >= 0.");
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
