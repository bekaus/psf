/*$Id: PeakShapeFunction.h 2496 2009-10-03 00:14:35Z mkirchner $*/

/*
 * PeakShapeFunction.h
 *
 * Copyright (c) 2008 Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
 *               2009 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
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

#ifndef __PEAKSHAPEFUNCTION_H__
#define __PEAKSHAPEFUNCTION_H__

#include <cmath>
#include <string>

#include <ms++/config.h>
#include <ms++/Error.h>
#include <ms++/Log.h>
#include <ms++/SparseSpectrum.h>

#include "psf/PeakParameter.h"
#include "psf/PeakShape.h"
#include "psf/PeakShapeFunctionTemplate.h"



/**
 * @page peakshapefunction Peakshape Functions
 *
 *
 *
 * @section introductionofpeakshapefunction What is a Peakshape Function?
 *
 * A mass spectrum usually consists of intensities drawn over the mass-over-charge values,
 * called 'mass channels'.
 * In theory, ions with a distinct mass-over-charge ratio should only contribute to one
 * mass channel producing a sharp peak in the mass spectrum. Due to imprecisions in the
 * measurement caused by random physical processes the actual peak is smeared out over a
 * few neighboring mass channels, resembling roughly a bell-shaped-curve in the intensity
 * domain.
 * Thus, a function describing the envelope of a single peak in a mass spectrum is called
 * a peakshape function.
 *
 *
 * @section typesofpeakshapefunctions Types of Peakshape Functions in ms++
 *
 * @subsection boxpsf Box
 * The box peakshape function is the rectangular function with a fixed width. It is implemented and used in the
 * transmogrifier unit test.
 *
 * This function is too crude an approximation and shouldn't be used in serious data analysis.
 *
 *
 * @subsection gaussianpsf Gaussian
 * The gaussian peakshape function is @f$ e^{-\frac{(\Delta\frac{m}{z})^{2}}{2\cdot\sigma^{2}}} @f$.
 * The only parameter is the full width at half maximum, which is connected to @f$ \sigma @f$ via
 * @f$ \mathrm{FWHM}=2\sqrt{2\ln2}\cdot\sigma @f$.
 *
 * This peakshape function is independent of the absolute value of the mass channel. In a real mass spectrum
 * the width of the peakshape is dependent on the mass channel (the resolution is usually decreasing
 * with higher masses). So, this function should only be used in a relatively small mass interval, where one
 * can safely neglect this change in the peak width.
 *
 * @see ms::GaussianPeakShapeFunction
 *
 *
 * @subsection orbipsf Orbitrap
 * The Orbitrap peak shape is gaussian. The full width at half maximum of the gaussian depends on the
 * absolute value of the mass channel and is calculated as follows: @f$ \mathrm{FWHM} = a\cdot\mathrm{mass}\sqrt{\mathrm{mass}} @f$.
 *
 * @see ms::OrbitrapPeakShapeFunction
 *
 *
 * @subsection orbiboxpsf Orbitrap Box
 * The Orbitrap box peakshape function is the rectangular function. The width of the rectangle depends on the
 * absolute value of the mass channel: @f$ \mathrm{width} = \mathrm{mass} \cdot \mathrm{sigma\_ppm}@f$.
 *
 * @f$ \mathrm{sigma\_ppm} @f$ is the only parameter and can be choosen as the resolution of an Orbitrap in parts-per-million.
 * Since for a real Orbitrap the peak width depends on the mass channel in a non-linear way and the peak shape is
 * far away from a box, this is a very crude approximation. Nevertheless, it can be used to speed up the calculation on good data sets
 * (high signal-to-noise ratio, few overlapping isotope patterns etc.) and can perform passable on such data.
 *
 * @see ms::OrbiBoxPeakShapeFunction
 *
 *
 * @subsection tofpsf Time-of-Flight
 * The time-of-flight peakshape function is gaussian. The full width at half maximum depends on the absolute value
 * of the mass channel as follows: @f$ \mathrm{FWHM} = a \cdot \sqrt{\mathrm{mass}} + b @f$.
 * @f$ a @f$ and @f$ b @f$ are parameters. They can be determined with the help of a peakshape function
 * estimator (by fitting the function on a given spectrum).
 * This formula can be deduced from the physics of a time-of-flight mass spectrometer.
 * 
 * @see ms::TOFPeakShapeFunction and ms::TOFPeakShapeFunctionEstimator
 *
 *
 *
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-06-09
 */

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

#endif /*__PEAKSHAPEFUNCTION_H__*/

