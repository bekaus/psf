/*$Id: PeakShapeFunctions-test.cpp 2617 2009-10-28 08:54:59Z bkausler $*/

/*
 * PeakShapeFunctions-test.cpp
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

#include <ms++/config.h>

#include <ms++/Error.h>
#include "ms++/Log.h"

#include "ms++/PeakShapeFunction.h"
#include "ms++/PeakShapeFunctions.h"

#include "testdata.h"

#include "unittest.hxx"

struct PsfTestSuite : vigra::test_suite {
    PsfTestSuite() : vigra::test_suite("PeakShapeFunctions") {
        add( testCase(&PsfTestSuite::testOrbitrapPeakShapeFunction) );
        add( testCase(&PsfTestSuite::testGaussianPeakShapeFunction) );
    }

    void testOrbitrapPeakShapeFunction() {
        ms::OrbitrapPeakShapeFunction orbi_psf;
        shouldEqual(orbi_psf.getType().toEnum(), ms::orbi);

        // construction
        ms::OrbitrapPeakShapeFunction psf1(0.1214);
        shouldEqual(psf1.getA(), 0.1214);

        // operator()
        orbi_psf.setA(0.0123);
        shouldEqualTolerance(orbi_psf(400., 402.), 0.998856, 0.000001); 
        shouldEqualTolerance(orbi_psf(400., 397.64), 0.998407, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 602.), 0.999661, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 597.64), 0.999528, 0.000001);

        // test half maximum at full width
        orbi_psf.setA(1.); // corresponds to const. fwhm of 2.0
        double fullMaximum = orbi_psf(1., 1.);
        double halfMaximum = orbi_psf(1., 1.5);
        shouldEqual(fullMaximum, 2. * halfMaximum);

        // auto calibration
        MSPP_LOG(ms::logINFO) << "Calibrating orbitrap peak shape function.";
        ms::SparseSpectrum spectrum(dirTestdata + "/PeakShapeFunctions/realistic_ms1.wsv");
        
        orbi_psf.setA(0); // reset 
        orbi_psf.calibrateFor(spectrum.begin(), spectrum.end());
        shouldEqualTolerance(orbi_psf.getA(), 1.19781e-06, 0.00001);       
    }

    void testGaussianPeakShapeFunction() {
        ms::GaussianPeakShapeFunction gaussian_psf;
        shouldEqual(gaussian_psf.getType().toEnum(), ms::gaussian);

        // construction
        ms::GaussianPeakShapeFunction psf(0.11442);
        shouldEqual(psf.getA(), 0.11442);

        // gaussian psf should be independent of the mass channel observed
        gaussian_psf.setA(3.);
        shouldEqualTolerance(gaussian_psf(400.,402.), 0.291632, 0.000001);
        shouldEqualTolerance(gaussian_psf(400., 397.64), 0.17982, 0.00001);
        shouldEqualTolerance(gaussian_psf(600.,602.), 0.291632, 0.000001);
        shouldEqualTolerance(gaussian_psf(600., 597.64), 0.17982, 0.00001);

        // test half maximum at full width
        gaussian_psf.setA(2.); // corresponds to const. fwhm of 2.0
        double fullMaximum = gaussian_psf(400., 400.);
        double halfMaximum = gaussian_psf(400., 401.);
        shouldEqual(fullMaximum, 2. * halfMaximum);

        // auto calibration
        MSPP_LOG(ms::logINFO) << "Calibrating gaussian peak shape function.";
        ms::SparseSpectrum spectrum(dirTestdata + "/PeakShapeFunctions/realistic_ms1.wsv");
        
        gaussian_psf.setA(0); // reset 
        gaussian_psf.calibrateFor(spectrum.begin(), spectrum.end());
        shouldEqualTolerance(gaussian_psf.getA(), 0.031325, 0.000001);
    }
};


int main()
{
    PsfTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


