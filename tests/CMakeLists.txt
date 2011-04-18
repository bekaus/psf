INCLUDE_DIRECTORIES(
    ${MSPP_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/testdata.h.in
        ${CMAKE_CURRENT_SOURCE_DIR}/testdata.h
        @ONLY IMMEDIATE
    )

#ADD_SUBDIRECTORY(testdata)

#### Sources
SET(SRCS_SPECTRALPEAK SpectralPeak-test.cpp)
SET(SRCS_SPECTRUMALGORITHM SpectrumAlgorithm-test.cpp)
SET(SRCS_PEAKSHAPEFUNCTIONTEMPLATE PeakShapeFunctionTemplate-test.cpp)
SET(SRCS_PEAKPARAMETER PeakParameter-test.cpp)
SET(SRCS_PEAKSHAPE PeakShape-test.cpp)
SET(SRCS_PEAKSHAPEFUNCTIONS PeakShapeFunctions-test.cpp)
SET(SRCS_PEAKSHAPEFUNCTION  PeakShapeFunction-test.cpp)

MACRO(ADD_MSPP_TEST name exe src)
    STRING(REGEX REPLACE "test_([^ ]+).*" "\\1" test "${exe}" )

    #build the test
    ADD_EXECUTABLE(${exe} ${src})
    #link the test
    TARGET_LINK_LIBRARIES(${exe} ms++)
    
    #add test to global list of unit test
    ADD_TEST(${name} ${exe})

    #Add target for the test
    STRING(REGEX REPLACE "test_([^ ]+).*" "unit_\\1" unittest_target "${exe}" )
    ADD_CUSTOM_TARGET(${unittest_target} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${exe})
ENDMACRO(ADD_MSPP_TEST name exe src)


#### Unit tests

#ADD_MSPP_TEST("Collection" test_collection ${SRCS_COLLECTION})

#ADD_MSPP_TEST("ErrorReporting" test_error ${SRCS_ERROR})

#ADD_MSPP_TEST("Logging" test_log ${SRCS_LOG})
ADD_MSPP_TEST("PeakParameter" test_peakparameter ${SRCS_PEAKPARAMETER})
ADD_MSPP_TEST("PeakShape" test_peakshape ${SRCS_PEAKSHAPE})
ADD_MSPP_TEST("PeakShapeFunction" test_peakshapefunction ${SRCS_PEAKSHAPEFUNCTION})
ADD_MSPP_TEST("PeakShapeFunctions" test_peakshapefunctions ${SRCS_PEAKSHAPEFUNCTIONS})
ADD_MSPP_TEST("PeakShapeFunctionTemplate" test_peakshapefunctiontemplate ${SRCS_PEAKSHAPEFUNCTIONTEMPLATE})
#ADD_MSPP_TEST("SparseSpectrum" test_sparsespectrum ${SRCS_SPARSESPECTRUM})
#ADD_MSPP_TEST("SparseSpectrumCollection" test_sparsespectrumcollection ${SRCS_SPARSESPECTRUMCOLLECTION})

ADD_MSPP_TEST("SpectralPeak" test_spectralpeak ${SRCS_SPECTRALPEAK})
ADD_MSPP_TEST("SpectrumAlgorithm" test_spectrumalgorithm ${SRCS_SPECTRUMALGORITHM})

