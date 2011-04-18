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

#include "ms++/PeakShapeFunction.h"
#include "psf/PeakShapeFunctions.h"


struct PsfTestSuite : vigra::test_suite {
    PsfTestSuite() : vigra::test_suite("PeakShapeFunction") {
        add( testCase(&PsfTestSuite::testPsfTypeEnum) );
        add( testCase(&PsfTestSuite::testPeakShapeFunctionType) );
        add( testCase(&PsfTestSuite::testGetType) );
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
};


int main()
{
    PsfTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


