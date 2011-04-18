/*$Id: PeakShapeFunctionTemplate-test.cpp 2306 2009-09-01 12:33:11Z bkausler $*/

/*
 * PeakShapeFunctionTemplate-test.cpp
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

#include <iostream>
#include <limits>

#include <ms++/Log.h>
#include <ms++/PeakParameter.h>
#include <ms++/PeakShape.h>
#include <ms++/PeakShapeFunction.h>
#include <ms++/PeakShapeFunctionTemplate.h>

#include "unittest.hxx"

#include "testdata.h"

using namespace ms;

struct PeakShapeFunctionTemplateTestSuite : vigra::test_suite {
    PeakShapeFunctionTemplateTestSuite() : vigra::test_suite("PeakShapeFunctionTemplate") {
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testOperator));
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testGetSupportThreshold));
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testClone));
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testGetType));
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testVirtuality));
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testSet_GetMinimalPeakHeightForCalibration));
        add( testCase(&PeakShapeFunctionTemplateTestSuite::testOrbiFwhmLinearSqrtPeakShape));
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

    void testClone() {
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::TofFwhm, ms::tof> gen_original;
        gen_original.setA(1.23456);
        gen_original.setB(6.54321);
        
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::TofFwhm, ms::tof>* cloned = gen_original.clone();
        shouldEqual(cloned->getA(), gen_original.getA());
        shouldEqual(cloned->getB(), gen_original.getB());

        delete cloned;
        cloned = NULL;   
    }

    void testGetType() {
        ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::TofFwhm, ms::tof> gen;
        shouldEqual(gen.getType().toEnum(), ms::tof);
    }

    void testVirtuality() {
        // we test the correct behaviour treating the Template as a PeakShapeFunction
        ms::PeakShapeFunction* psf = new ms::PeakShapeFunctionTemplate<ms::GaussianPeakShape, ms::TofFwhm, ms::tof>;
        shouldEqual(psf->getType().toEnum(), ms::tof);
        
        delete psf;
        psf = NULL;
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
    PeakShapeFunctionTemplateTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


