/* $Id: PeakShapeFunctions.h 2496 2009-10-03 00:14:35Z mkirchner $ */

/*
 * PeakShapeFunctions.h
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

#ifndef __PEAKSHAPEFUNCTIONS_H__
#define __PEAKSHAPEFUNCTIONS_H__

#include <ms++/PeakParameter.h>
#include <ms++/PeakShape.h>
#include <ms++/PeakShapeFunction.h>
#include <ms++/PeakShapeFunctionTemplate.h>

namespace ms
{

/**
* A peak shape function as it occurs in Orbitrap mass spectra.
*
* The Orbitrap peak shape function is parameterized via a linear sqrt model, which goes through the origin: @f$ f(x) = a\cdot x\sqrt{x}@f$ .
* You can set the model parameter a via the corresponding getter/setter method in the
* PeakShapeFunctionTemplate interface. Furthermore, you may autocalibrate the parameter calling the calibrateFor() method.
*
* This function is robust concerning autocalibration, because the autocalibration cannot set the parameter a so, that the function becomes invalid in some
* mz ranges. That's why we choose a PeakParameter, which is constrained in the origin.
* @see ms::RobustOrbitrapPeakShapeFunction 
*/
typedef PeakShapeFunctionTemplate<GaussianPeakShape, OrbitrapWithOriginFwhm, orbi> OrbitrapPeakShapeFunction;

/**
* A peak shape function as it occurs in centroided Orbitrap mass spectra.
*
* The function is similar to the @see ms::OrbitrapPeakShapeFunction, the only difference being that 
* the window shape is a box. Support threshold calculation etc are identical to the @see OrbitrapPeakShapeFunction.
* The Orbitrap peak shape function is parameterized via a linear sqrt model, which goes through the origin: @f$ f(x) = a\cdot x\sqrt{x}@f$ .
* You can set the model parameter a via the corresponding getter/setter method in the
* PeakShapeFunctionTemplate interface. Furthermore, you may autocalibrate the parameter calling the calibrateFor() method.
*
* This function is robust concerning autocalibration, because the autocalibration cannot set the parameter a so, that the function becomes invalid in some
* mz ranges. That's why we choose a PeakParameter, which is constrained in the origin.
* @see ms::RobustOrbitrapPeakShapeFunction 
*/
typedef PeakShapeFunctionTemplate<BoxPeakShape, OrbitrapWithOriginFwhm, orbi> OrbitrapBoxPeakShapeFunction;

/**
* A peak shape function with a gaussian shape static everywhere in a mass spectrum.
*
* You can set the full width at half maximum of the gaussian via the 'a' getter and setter 
* in the PeakShapeFunctionTemplate
* interface (the 'b' functions are not supported and produce compile time errors if used).
* Furthermore, you may autocalibrate the parameter calling the calibrateFor() method.
*/
typedef PeakShapeFunctionTemplate<GaussianPeakShape, ConstantFwhm, gaussian> GaussianPeakShapeFunction;

} /* namespace ms */

#endif /*__PEAKSHAPEFUNCTIONS_H__*/
