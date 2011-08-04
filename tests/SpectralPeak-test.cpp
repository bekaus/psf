/*$Id: SpectralPeak-test.cpp 2620 2009-10-28 09:42:57Z bkausler $*/

/*
 * SpectralPeak-test.cpp
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

#include <ms++/Error.h>
#include <psf/SpectralPeak.h>
#include <psf/Spectrum.h>

#include "unittest.hxx"

using namespace ms;

struct SpectralPeakTestSuite : vigra::test_suite {
    SpectralPeakTestSuite() : vigra::test_suite("SpectralPeak") {
        add( testCase(&SpectralPeakTestSuite::testHeight));
        add( testCase(&SpectralPeakTestSuite::testLowness));
        add( testCase(&SpectralPeakTestSuite::testFullWidthAtFractionOfMaximum));
    }

    void testHeight() {
	IntensityExtractor get_int;
        // A peak with height 3.1
        Spectrum s1;
        s1.push_back(SpectrumElement(1.1, 1.1));
        s1.push_back(SpectrumElement(1.2, 1.9));
        s1.push_back(SpectrumElement(1.4, 3.1));
        s1.push_back(SpectrumElement(1.5, 2.2));
        s1.push_back(SpectrumElement(1.69, 1.14));
        s1.push_back(SpectrumElement(1.76, 0.98));
	shouldEqual(SpectralPeak::height(get_int, s1.begin(), --(s1.end())), 3.1);
        
        // A empty spectrum
        Spectrum s2;
        bool thrown = false;
        try {
            SpectralPeak::height(get_int, s2.begin(), --(s2.end()));
        }
        catch(const PreconditionViolation& e) {
			MSPP_UNUSED(e);
            thrown = true;
        }
        should(thrown);
    }

    void testLowness() {
	IntensityExtractor get_int;
        // A quite normal spectral peak.
        //
        // The maximum abundance is 3.1
        // The lowest abundance on the 'left' is 1.1 and on the 'right' 0.98.
        // So, the lowness is (1 - 1.1/3.1) = 0.64516129.
        Spectrum s1;
        s1.push_back(SpectrumElement(1.1, 1.1));
        s1.push_back(SpectrumElement(1.2, 1.9));
        s1.push_back(SpectrumElement(1.4, 3.1));
        s1.push_back(SpectrumElement(1.5, 2.2));
        s1.push_back(SpectrumElement(1.69, 1.14));
        s1.push_back(SpectrumElement(1.76, 0.98));

        shouldEqual(SpectralPeak::lowness(get_int, s1.begin(), --(s1.end())), 1-(1.1/3.1));      

        // A peak with only one 'flank'.
        //
        // Lowness of a one flanked peak is 0.0.
        Spectrum s5;
        s5.push_back(SpectrumElement(1.1, 1.1));
        s5.push_back(SpectrumElement(1.2, 1.9));
        s5.push_back(SpectrumElement(1.4, 3.1));
        s5.push_back(SpectrumElement(1.5, 5.2));

        shouldEqual(SpectralPeak::lowness(get_int, s5.begin(), --(s5.end())), 0.0);     

        // An equiabundant sequence.
        //
        // The lowness of an equiabundant sequence is 0.0.
        Spectrum s2;
        s2.push_back(SpectrumElement(1.1, 1.1));
        s2.push_back(SpectrumElement(1.2, 1.1));
        s2.push_back(SpectrumElement(1.4, 1.1));
        s2.push_back(SpectrumElement(1.5, 1.1));

        shouldEqual(SpectralPeak::lowness(get_int, s2.begin(), --(s2.end())), 0.0);  

        // Quite unrealistic peak, with zero abundance elements
        //
        // Lowness is than a straight 1.0.
        Spectrum s6;
        s6.push_back(SpectrumElement(1.1, 0.1));
        s6.push_back(SpectrumElement(1.2, 0.0));
        s6.push_back(SpectrumElement(1.4, 1.1));
        s6.push_back(SpectrumElement(1.5, 1.2));
        s6.push_back(SpectrumElement(1.7, 0.0));
        s6.push_back(SpectrumElement(1.9, 1.1));
        s6.push_back(SpectrumElement(2.12, 0.9));

        shouldEqual(SpectralPeak::lowness(get_int, s6.begin(), --(s6.end())), 1.0);

        // Only one element.
        //
        // The lowness of one element is 0.0.
        Spectrum s3;
        s3.push_back(SpectrumElement(123.32, 89.1));

        shouldEqual(SpectralPeak::lowness(get_int, s3.begin(), --(s3.end())), 0.0);
    }

    void testFullWidthAtFractionOfMaximum() {
        // A 'normal' peak.
        // Note the abundance twist in the last two elements.
        //
        // Fraction | Full width
        // 0.7      | 0.257459
        // 0.5      | 0.397029
        // 0.3      | lowness to small -> not defined     
        SparseSpectrum s1;
        s1.push_back(SparseSpectrum::Element(0.4, 0.12));
        s1.push_back(SparseSpectrum::Element(1.1, 1.1));
        s1.push_back(SparseSpectrum::Element(1.2, 1.9));
        s1.push_back(SparseSpectrum::Element(1.4, 3.1));
        s1.push_back(SparseSpectrum::Element(1.5, 2.2));
        s1.push_back(SparseSpectrum::Element(1.6, 0.98));
        s1.push_back(SparseSpectrum::Element(1.69, 1.14));

        shouldEqualTolerance(SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), 0.7), 0.257459, 0.000001);
        shouldEqualTolerance(SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), 0.5), 0.397029, 0.000001);

        // not defined: full width at fraction of 0.3
        bool thrown = false;
        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), 0.3);
        }
        catch(const Starvation& e) {
			MSPP_UNUSED(e);
            thrown = true;        
        }
        shouldMsg(thrown, "Starvation wasn't thrown despite of invalid input data.");
        thrown = false;

        // illegal fraction
        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), 1.1);
        }
        catch(const PreconditionViolation& e) {
			MSPP_UNUSED(e);
            thrown = true;        
        }
        shouldMsg(thrown, "PreconditionViolation wasn't thrown despite of invalid fraction parameter.");
        thrown = false;

        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), -0.3);
        }
        catch(const PreconditionViolation& e) {
			MSPP_UNUSED(e);
            thrown = true;        
        }
        shouldMsg(thrown, "PreconditionViolation wasn't thrown despite of invalid fraction parameter.");
        thrown = false;

        // legal border fraction values
        // no PreconditionViolation should be thrown here
        SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), 1.0);
        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(s1.begin(), --(s1.end()), 0.0);
        }
        catch(const Starvation& e) {
            MSPP_UNUSED(e);
        }

        // Special peak with elements exactly on target abundance
        SparseSpectrum s_onTarget;
        s_onTarget.push_back(SparseSpectrum::Element(3, 7));
        s_onTarget.push_back(SparseSpectrum::Element(4, 10));
        s_onTarget.push_back(SparseSpectrum::Element(5, 7));

        shouldEqualTolerance(SpectralPeak::fullWidthAtFractionOfMaximum(s_onTarget.begin(), --(s_onTarget.end()), 0.71), 2., 0.1);
    }
};

int main()
{
    SpectralPeakTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


