/*$Id: PeakShape-test.cpp 2620 2009-10-28 09:42:57Z bkausler $*/

/*
 * PeakShape-test.cpp
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
#include <iostream>

#include <ms++/Error.h>
#include "psf/PeakShape.h"

#include "unittest.hxx"


struct peakshapeTestSuite : vigra::test_suite {
    peakshapeTestSuite() : vigra::test_suite("peakshape") {
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeConstruction));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeGetterSetter));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeSigmaFwhmConversion));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeAt));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeGetSupportThreshold));
    }

    void testGaussianPeakShapeConstruction() {
        // default constructor
        psf::GaussianPeakShape gps;
        shouldEqual(gps.getSigma(), 0.1);

        // constructor init sigma
        psf::GaussianPeakShape gps_sigma(0.79);
        shouldEqual(gps_sigma.getSigma(), 0.79);

        // illegal sigma
        bool thrown = false;
        try {
            psf::GaussianPeakShape(0);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            psf::GaussianPeakShape(-0.34);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;
    }

    void testGaussianPeakShapeGetterSetter() {
        psf::GaussianPeakShape gps;
        bool thrown = false;

        // sigma
        gps.setSigma(0.5);
        shouldEqual(gps.getSigma(), 0.5);
                
        gps.setSigma(0.7);
        shouldEqual(gps.getSigma(), 0.7);

        try {
            gps.setSigma(0.);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            gps.setSigma(-1.7);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // fwhm
        gps.setFwhm(0.5);
        shouldEqual(gps.getFwhm(), 0.5);

        gps.setFwhm(0.7);
        shouldEqual(gps.getFwhm(), 0.7);

        try {
            gps.setFwhm(0.);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            gps.setFwhm(-1.7);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // sigmaFactorForSupportThreshold
        gps.setSigmaFactorForSupportThreshold(0.5);
        shouldEqual(gps.getSigmaFactorForSupportThreshold(), 0.5);

        gps.setSigmaFactorForSupportThreshold(0.7);
        shouldEqual(gps.getSigmaFactorForSupportThreshold(), 0.7);

        try {
            gps.setSigmaFactorForSupportThreshold(0.);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            gps.setSigmaFactorForSupportThreshold(-1.7);
        } catch (const psf::PreconditionViolation& e){
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;
    }

    void testGaussianPeakShapeSigmaFwhmConversion() {
        psf::GaussianPeakShape gps;

        // conversion factor
        double conversionFactor = gps.sigmaToFwhmConversionFactor();
        shouldEqualTolerance(conversionFactor, 2.35482, 0.00001);

        // sigma / fwhm conversion
        gps.setSigma(0.5);
        shouldEqual(gps.getFwhm(), conversionFactor * 0.5);

        gps.setFwhm(0.5);
        shouldEqual(gps.getSigma(), 0.5 / conversionFactor);
    }

    void testGaussianPeakShapeAt() {
        struct Gauss {
            static double at(const double x, const double sigma) {
                return std::exp(-(x * x) / (2 * sigma * sigma));
            }
        };

        psf::GaussianPeakShape gps;

        // has to be 1 at xCoordinate=0 
        shouldEqual(gps.at(0), 1);

        // test some values
        gps.setSigma(0.5);
        shouldEqual(gps.at(0.1), Gauss::at(0.1, 0.5));
        shouldEqual(gps.at(3.5), Gauss::at(3.5, 0.5));
        shouldEqual(gps.at(-0.34), Gauss::at(-0.34, 0.5));
        shouldEqual(gps.at(-2.73), Gauss::at(-2.73, 0.5));

        gps.setSigma(0.9);
        shouldEqual(gps.at(0.1), Gauss::at(0.1, 0.9));
        shouldEqual(gps.at(3.5), Gauss::at(3.5, 0.9));
        shouldEqual(gps.at(-0.34), Gauss::at(-0.34, 0.9));
        shouldEqual(gps.at(-2.73), Gauss::at(-2.73, 0.9));
    }

    void testGaussianPeakShapeGetSupportThreshold() {
        psf::GaussianPeakShape gps;

        shouldEqual(gps.getSigmaFactorForSupportThreshold(), 3.0);

        gps.setSigma(1.5);
        shouldEqual(gps.getSupportThreshold(), 4.5);

        gps.setSigma(0.7);
        shouldEqual(gps.getSupportThreshold(), 0.7 * 3.0);
    }
};

int main()
{
    peakshapeTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


