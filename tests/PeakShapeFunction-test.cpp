/*$Id: PeakShapeFunction-test.cpp 2253 2009-08-19 11:43:25Z bkausler $*/

/*
 * PeakShapeFunction-test.cpp
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

#include "unittest.hxx"

#include "ms++/Log.h"

#include "psf/PeakShapeFunction.h"
#include <psf/Spectrum.h>
#include "testdata.h"

struct PsfTestSuite : vigra::test_suite {
    PsfTestSuite() : vigra::test_suite("PeakShapeFunction") {
        add( testCase(&PsfTestSuite::testPsfTypeEnum) );
        add( testCase(&PsfTestSuite::testPeakShapeFunctionType) );
        add( testCase(&PsfTestSuite::testGetType) );
        add( testCase(&PsfTestSuite::testOrbitrapPeakShapeFunction) );
        add( testCase(&PsfTestSuite::testGaussianPeakShapeFunction) );
        add( testCase(&PsfTestSuite::testOperator));
        add( testCase(&PsfTestSuite::testGetSupportThreshold));
        add( testCase(&PsfTestSuite::testSet_GetMinimalPeakHeightForCalibration));
        add( testCase(&PsfTestSuite::testOrbiFwhmLinearSqrtPeakShape));
    }


    void testPsfTypeEnum() {
        // there should be six different enum values: box, gaussian, orbi, orbiBox, orbiConst, tof
        MSPP_LOG(ms::logINFO) << "Testing the PeakShapeFunctionType enum.";
        ms::PeakShapeFunctionTypes t;
        t = ms::box;
        t = ms::gaussian;
        t = ms::orbi;
        t = ms::orbiBox;
        t = ms::tof;
    }


    void testPeakShapeFunctionType() {
        MSPP_LOG(ms::logINFO) << "Testing class PeakShapeFunctionType";

        ms::PeakShapeFunctionType psfBox = ms::box;
        shouldEqual(psfBox.toEnum(), ms::box);
        should(psfBox.toString() == "box");

        ms::PeakShapeFunctionType psfGaussian = ms::gaussian;
        shouldEqual(psfGaussian.toEnum(), ms::gaussian);
        should(psfGaussian.toString() == "gaussian");

        ms::PeakShapeFunctionType psfOrbi = ms::orbi;
        shouldEqual(psfOrbi.toEnum(), ms::orbi);
        should(psfOrbi.toString() == "orbi");

        ms::PeakShapeFunctionType psfOrbiBox = ms::orbiBox;
        shouldEqual(psfOrbiBox.toEnum(), ms::orbiBox);
        should(psfOrbiBox.toString() == "orbiBox");

        ms::PeakShapeFunctionType psfTof = ms::tof;
        shouldEqual(psfTof.toEnum(), ms::tof);
        should(psfTof.toString() == "time-of-flight");

        // test illegal enum (choose the integer high enough...)
        ms::PeakShapeFunctionType psfIllegal = static_cast<ms::PeakShapeFunctionTypes>(200);
        shouldEqual(psfIllegal.toEnum(), static_cast<ms::PeakShapeFunctionTypes>(200));
        should(psfIllegal.toString() == "unknown");

    }
    
    // We want to test this function for every implementation of the abstract
    // 'PeakShapeFunction' interface.
    // The 'box' type is only used for unit testing and not included in the core
    // library. So, we don't test it.
    void testGetType() {
        // Testing getType(). We use sound values for the constructors.
        MSPP_LOG(ms::logINFO) << "Testing the getType() functions.";
        shouldEqual(ms::GaussianPeakShapeFunction().getType().toEnum(), ms::gaussian);
        shouldEqual(ms::OrbitrapPeakShapeFunction().getType().toEnum(), ms::orbi);
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
	ms::Spectrum spectrum;
	loadSpectrumElements(spectrum,dirTestdata + "/PeakShapeFunctions/realistic_ms1.wsv");
	ms::MzExtractor get_mz;
	ms::IntensityExtractor get_int;
        
        orbi_psf.setA(0); // reset 
        orbi_psf.calibrateFor(get_mz, get_int, spectrum.begin(), spectrum.end());
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
	ms::Spectrum spectrum;
	loadSpectrumElements(spectrum,dirTestdata + "/PeakShapeFunctions/realistic_ms1.wsv");
	ms::MzExtractor get_mz;
	ms::IntensityExtractor get_int;
        
        gaussian_psf.setA(0); // reset 
        gaussian_psf.calibrateFor(get_mz, get_int, spectrum.begin(), spectrum.end());
        shouldEqualTolerance(gaussian_psf.getA(), 0.031325, 0.000001);
    }

    void testOperator() {
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::TofFwhm, ms::tof> gen;
        ms::TofFwhm fwhm;
        ms::GaussianPeakShape ps;

        gen.setA(0.43);
        gen.setB(0.76);
        fwhm.setA(0.43);
        fwhm.setB(0.76);

        ps.setFwhm(fwhm.at(400.));

        // ensure to be inside the threshold
        should(gen.getSupportThreshold(400.) > 5.0);

        shouldEqual(gen(400., 404.5), ps.at(404.5 - 400.));
        shouldEqual(gen(400., 397.2), ps.at(400. - 397.2));
        shouldEqual(gen(400., 400.), ps.at(0.));

        // now, test the threshold behaviour

        double threshold = gen.getSupportThreshold(400.);
        // Should be not zero at the threshold
        should(gen(400., 400. + threshold) > 0.);
        should(gen(400., 400. - threshold) > 0.);

        // zero after the threshold
        double delta = std::numeric_limits<double>().epsilon() + std::numeric_limits<double>().round_error();
        shouldEqual(gen(400., 400. + threshold + delta), 0.0);
        shouldEqual(gen(400., 400. - (threshold + delta)), 0.0);       
    }

    void testGetSupportThreshold() {
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::TofFwhm, ms::tof> gen;
        ms::TofFwhm fwhm;
        ms::GaussianPeakShape ps;

        gen.setA(0.43);
        gen.setB(0.76);
        fwhm.setA(0.43);
        fwhm.setB(0.76);

        ps.setFwhm(fwhm.at(400.));
        
        double threshold = ps.getSupportThreshold();
        shouldEqual(gen.getSupportThreshold(400.), threshold);
    }

    void testSet_GetMinimalPeakHeightForCalibration() {
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::OrbitrapFwhm, ms::orbi> psf;
       
        psf.setMinimalPeakHeightForCalibration(4.2);
        shouldEqual(psf.getMinimalPeakHeightForCalibration(), 4.2);

        psf.setMinimalPeakHeightForCalibration(0);
        shouldEqual(psf.getMinimalPeakHeightForCalibration(), 0);

        psf.setMinimalPeakHeightForCalibration(-0.87);
        shouldEqual(psf.getMinimalPeakHeightForCalibration(), -0.87);
    }

    void testOrbiFwhmLinearSqrtPeakShape() {
        // Note: This was the test for the first version of the orbi psf, but doesn't apply anymore.
        // Nevertheless, it is still a valid test for its combination of template parameters. Therefore,
        // we 'retired' the test here.

        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::OrbitrapFwhm, ms::orbi>  orbi_psf;
        shouldEqual(orbi_psf.getType().toEnum(), ms::orbi);

        // construction
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::OrbitrapFwhm, ms::orbi> psf1(0.1214);
        shouldEqual(psf1.getA(), 0.1214);
        shouldEqual(psf1.getB(), orbi_psf.getB());

        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::OrbitrapFwhm, ms::orbi> psf2(0.1214, 1.20392);
        shouldEqual(psf2.getA(), 0.1214);
        shouldEqual(psf2.getB(), 1.20392);

        // operator()
        orbi_psf.setA(0.0123);
        orbi_psf.setB(0.0234);
        shouldEqualTolerance(orbi_psf(400., 402.), 0.998856, 0.000001); 
        shouldEqualTolerance(orbi_psf(400., 397.64), 0.998407, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 602.), 0.999661, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 597.64), 0.999528, 0.000001);

        // test half maximum at full width
        orbi_psf.setA(1.); // corresponds to const. fwhm of 2.0
        orbi_psf.setB(0.);
        double fullMaximum = orbi_psf(1., 1.);
        double halfMaximum = orbi_psf(1., 1.5);
        shouldEqual(fullMaximum, 2. * halfMaximum); 
    }

};


int main()
{
    PsfTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


