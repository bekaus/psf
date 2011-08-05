/*$Id: PeakParameter-test.cpp 2620 2009-10-28 09:42:57Z bkausler $*/

/*
 * PeakParameter-test.cpp
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

#include <ms++/config.h>
#include <ms++/Error.h>
#include <psf/PeakParameter.h>
#include <psf/Spectrum.h>

#include "testdata.h"

#include "unittest.hxx"

struct PeakParameterTestSuite : vigra::test_suite {
    PeakParameterTestSuite() : vigra::test_suite("PeakParameter") {
        add( testCase(&PeakParameterTestSuite::testSet_GetMinimalPeakHeightToLearnFrom));
        add( testCase(&PeakParameterTestSuite::testOrbitrapFwhm));
        add( testCase(&PeakParameterTestSuite::testFtIcrFwhm));
        add( testCase(&PeakParameterTestSuite::testTofFwhm));
        add( testCase(&PeakParameterTestSuite::testConstantFwhm));
        add( testCase(&PeakParameterTestSuite::testConstantFwhmLearnFrom));
        add( testCase(&PeakParameterTestSuite::testOrbitrapFwhmLearnFrom));
        add( testCase(&PeakParameterTestSuite::testFtIcrFwhmLearnFrom));
        add( testCase(&PeakParameterTestSuite::testTofFwhmLearnFrom));
    }

    void testSet_GetMinimalPeakHeightToLearnFrom() {
        using namespace psf;
        PeakParameterFwhm<ConstantModel> fwhm;

        // should be zero by default after construction
        shouldEqual(fwhm.getMinimalPeakHeightToLearnFrom(), 0);
        
        fwhm.setMinimalPeakHeightToLearnFrom(0.92);
        shouldEqual(fwhm.getMinimalPeakHeightToLearnFrom(), 0.92);

        fwhm.setMinimalPeakHeightToLearnFrom(0);
        shouldEqual(fwhm.getMinimalPeakHeightToLearnFrom(), 0);

        fwhm.setMinimalPeakHeightToLearnFrom(-1.7);
        shouldEqual(fwhm.getMinimalPeakHeightToLearnFrom(), -1.7);
    }

    void testOrbitrapFwhm() {
        psf::OrbitrapFwhm fwhm;
    
        // number of parameters
        shouldEqual(fwhm.numberOfParameters(), unsigned(2));

        // setter / getter
        fwhm.setA(234.3);
        shouldEqual(fwhm.getA(), 234.3);
        fwhm.setA(-234.321);
        shouldEqual(fwhm.getA(), -234.321);
        fwhm.setA(0);
        shouldEqual(fwhm.getA(), 0);

        fwhm.setB(234.3);
        shouldEqual(fwhm.getB(), 234.3);
        fwhm.setB(-234.321);
        shouldEqual(fwhm.getB(), -234.321);
        fwhm.setB(0);
        shouldEqual(fwhm.getB(), 0);

        fwhm.setParameter(0, 9437);
        shouldEqual(fwhm.getParameter(0), 9437);
        fwhm.setParameter(0, -9437.1);
        shouldEqual(fwhm.getParameter(0), -9437.1);
        fwhm.setParameter(0, 0);
        shouldEqual(fwhm.getParameter(0), 0);

        fwhm.setParameter(1, 9437.1);
        shouldEqual(fwhm.getParameter(1), 9437.1);
        fwhm.setParameter(1, -9437.1);
        shouldEqual(fwhm.getParameter(1), -9437.1);
        fwhm.setParameter(1, 0);
        shouldEqual(fwhm.getParameter(1), 0);

        bool thrown = false;
        try {
            fwhm.setParameter(2,0);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        thrown = false;
        try {
            fwhm.getParameter(2);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // at()      
        fwhm.setA(0.43);
        fwhm.setB(0.76);
        shouldEqualTolerance(fwhm.at(400), 3440.76, 1e-2);

        // no masses <= 0
        try {
            fwhm.at(-123.2);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        try {
            fwhm.at(0);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        // negative fwhm
        fwhm.setA(-0.1);
        fwhm.setB(0.1);
        try {
            fwhm.at(400);
        } catch (const psf::PostconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;      
    }


    void testFtIcrFwhm() {
        psf::FtIcrFwhm fwhm;

        // number of parameters
        shouldEqual(fwhm.numberOfParameters(), unsigned(2));

        // setter / getter
        fwhm.setA(234.3);
        shouldEqual(fwhm.getA(), 234.3);
        fwhm.setA(-234.321);
        shouldEqual(fwhm.getA(), -234.321);
        fwhm.setA(0);
        shouldEqual(fwhm.getA(), 0);

        fwhm.setB(234.3);
        shouldEqual(fwhm.getB(), 234.3);
        fwhm.setB(-234.321);
        shouldEqual(fwhm.getB(), -234.321);
        fwhm.setB(0);
        shouldEqual(fwhm.getB(), 0);

        fwhm.setParameter(0, 9437.1);
        shouldEqual(fwhm.getParameter(0), 9437.1);
        fwhm.setParameter(0, -9437.1);
        shouldEqual(fwhm.getParameter(0), -9437.1);
        fwhm.setParameter(0, 0);
        shouldEqual(fwhm.getParameter(0), 0);

        fwhm.setParameter(1, 9437.1);
        shouldEqual(fwhm.getParameter(1), 9437.1);
        fwhm.setParameter(1, -9437.1);
        shouldEqual(fwhm.getParameter(1), -9437.1);
        fwhm.setParameter(1, 0);
        shouldEqual(fwhm.getParameter(1), 0);

        bool thrown = false;
        try {
            fwhm.setParameter(2,0);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        thrown = false;
        try {
            fwhm.getParameter(2);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // at()      
        fwhm.setA(0.43);
        fwhm.setB(0.76);
        shouldEqual(fwhm.at(400), 68800.76);

        // no masses <= 0
        thrown = false;
        try {
            fwhm.at(-123.2);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        try {
            fwhm.at(0);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        // negative fwhm
        fwhm.setA(-0.1);
        fwhm.setB(0.1);
        try {
            fwhm.at(400);
        } catch (const psf::PostconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;
    }

    
    void testTofFwhm() {
        psf::TofFwhm fwhm;

        // number of parameters
        shouldEqual(fwhm.numberOfParameters(), unsigned(2));

        // setter / getter
        fwhm.setA(234.3);
        shouldEqual(fwhm.getA(), 234.3);
        fwhm.setA(-234.321);
        shouldEqual(fwhm.getA(), -234.321);
        fwhm.setA(0);
        shouldEqual(fwhm.getA(), 0);

        fwhm.setB(234.3);
        shouldEqual(fwhm.getB(), 234.3);
        fwhm.setB(-234.321);
        shouldEqual(fwhm.getB(), -234.321);
        fwhm.setB(0);
        shouldEqual(fwhm.getB(), 0);

        fwhm.setParameter(0, 9437.1);
        shouldEqual(fwhm.getParameter(0), 9437.1);
        fwhm.setParameter(0, -9437.1);
        shouldEqual(fwhm.getParameter(0), -9437.1);
        fwhm.setParameter(0, 0);
        shouldEqual(fwhm.getParameter(0), 0);

        fwhm.setParameter(1, 9437.1);
        shouldEqual(fwhm.getParameter(1), 9437.1);
        fwhm.setParameter(1, -9437.1);
        shouldEqual(fwhm.getParameter(1), -9437.1);
        fwhm.setParameter(1, 0);
        shouldEqual(fwhm.getParameter(1), 0);

        bool thrown = false;
        try {
            fwhm.setParameter(2,0);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        thrown = false;
        try {
            fwhm.getParameter(2);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // at()      
        fwhm.setA(0.43);
        fwhm.setB(0.76);
        shouldEqual(fwhm.at(400), 9.36);

        // no masses <= 0
        thrown = false;
        try {
            fwhm.at(-123.2);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        try {
            fwhm.at(0);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        // negative fwhm
        fwhm.setA(-0.1);
        fwhm.setB(0.1);
        try {
            fwhm.at(400);
        } catch (const psf::PostconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;
    }

    void testConstantFwhm() {
        psf::ConstantFwhm fwhm;

        // number of parameters
        shouldEqual(fwhm.numberOfParameters(), unsigned(1));

        // setter / getter
        fwhm.setA(234.3);
        shouldEqual(fwhm.getA(), 234.3);
        fwhm.setA(-234.321);
        shouldEqual(fwhm.getA(), -234.321);
        fwhm.setA(0);
        shouldEqual(fwhm.getA(), 0);

        fwhm.setParameter(0, 9437.1);
        shouldEqual(fwhm.getParameter(0), 9437.1);
        fwhm.setParameter(0, -9437.1);
        shouldEqual(fwhm.getParameter(0), -9437.1);
        fwhm.setParameter(0, 0);
        shouldEqual(fwhm.getParameter(0), 0);

        bool thrown = false;
        try {
            fwhm.setParameter(1,0);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        thrown = false;
        try {
            fwhm.getParameter(1);
        } catch(const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // at()
        fwhm.setA(0.43);
        shouldEqual(fwhm.at(100), fwhm.getA());
        shouldEqual(fwhm.at(400), fwhm.getA());
        
        /* dirty tests */        
        thrown = false;
        
        // negative and zero masses
        fwhm.setA(0.1);
        try {
            fwhm.at(-123.2);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;

        try {
            fwhm.at(0);
        } catch (const psf::PreconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;
        
        // negative fwhm
        fwhm.setA(-0.1);
        try {
            fwhm.at(400);
        } catch (const psf::PostconditionViolation& e) {
			PSF_UNUSED(e);
            thrown = true;        
        }
        should(thrown);
        thrown = false;
    }

    void testConstantFwhmLearnFrom() {
        using namespace psf;
        using namespace std;
	MzExtractor get_mz;
	IntensityExtractor get_int;
        PSF_LOG(logINFO) << "Testing ConstantFwhm learnFrom().";
        ConstantFwhm fwhm;
	Spectrum spectrum;
        loadSpectrumElements(spectrum, dirTestdata + "/PeakParameter/realistic_ms1.wsv");
        
        fwhm.setA(0); // reset A
        fwhm.learnFrom(get_mz, get_int, spectrum.begin(), spectrum.end());
        shouldEqualTolerance(fwhm.getA(), 0.031325, 0.000001);
    }

    void testOrbitrapFwhmLearnFrom() {
        using namespace psf;
        using namespace std;
	MzExtractor get_mz;
	IntensityExtractor get_int;
        PSF_LOG(logINFO) << "Testing OrbitrapFwhm learnFrom().";
        OrbitrapFwhm fwhm;
	Spectrum spectrum;
        loadSpectrumElements(spectrum, dirTestdata + "/shared_data/orbi_ms1.wsv");
        
        fwhm.setA(0); // reset 
        fwhm.setB(0);
        fwhm.learnFrom(get_mz, get_int, spectrum.begin(), spectrum.end());
        shouldEqualTolerance(fwhm.getA(), 9.40679e-06, 0.00001);
        shouldEqualTolerance(fwhm.getB(), 0., 0.);
    }

    void testFtIcrFwhmLearnFrom() {
        using namespace psf;
        using namespace std;
	MzExtractor get_mz;
	IntensityExtractor get_int;
        PSF_LOG(logINFO) << "Testing FtIcrFwhm learnFrom()."; 
        FtIcrFwhm fwhm;
	Spectrum spectrum;
        loadSpectrumElements(spectrum, dirTestdata + "/PeakParameter/realistic_ms1.wsv");

        
        fwhm.setA(0); // reset 
        fwhm.setB(0);
        fwhm.learnFrom(get_mz, get_int, spectrum.begin(), spectrum.end());
        PSF_LOG(logINFO) << "Learned FtIcrFwhm Parameter A: " << fwhm.getA();
        PSF_LOG(logINFO) << "Learned FtIcrFwhm Parameter B: " << fwhm.getB();
        shouldEqualTolerance(fwhm.getA(), 0., 0.00001);
        shouldEqualTolerance(fwhm.getB(), 0.031325, 0.00001);
    }

    void testTofFwhmLearnFrom() {
        using namespace psf;
        using namespace std;
	MzExtractor get_mz;
	IntensityExtractor get_int;
        PSF_LOG(logINFO) << "Testing TofFwhm learnFrom().";
        TofFwhm fwhm;
	Spectrum spectrum;
        loadSpectrumElements(spectrum, dirTestdata + "/PeakParameter/realistic_ms1.wsv");
        
        fwhm.setA(0); // reset 
        fwhm.setB(0);
        fwhm.learnFrom(get_mz, get_int, spectrum.begin(), spectrum.end());
        shouldEqualTolerance(fwhm.getA(), 0., 0.00001);
        shouldEqualTolerance(fwhm.getB(), 0.031325, 0.0001);
    }
};

int main()
{
    PeakParameterTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


