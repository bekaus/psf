#include <string>

#include "psf/PeakShapeFunction.h"

using namespace psf;

PeakShapeFunctionTypes PeakShapeFunctionType::toEnum() {
    return type_;
}

std::string PeakShapeFunctionType::toString() {
    switch(type_) {
        case box:
            return "box";
            break;

        case gaussian:
            return "gaussian";
            break;

        case orbi:
            return "orbi";
            break;

        case orbiBox:
            return "orbiBox";
            break;

        case tof:
            return "time-of-flight";
            break;

        default:
            return "unknown";
            break;
    }
}
