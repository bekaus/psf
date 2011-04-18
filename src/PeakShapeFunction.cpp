/*$Id: PeakShapeFunction.cpp 2670 2009-11-09 15:09:59Z bkausler $*/

/*
 * PeakShapeFunction.cpp
 *
 * Copyright (c) 2009 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * Copyright (c) 2008 Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
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

#include <string>

#include "psf/PeakShapeFunction.h"

using namespace ms;

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
