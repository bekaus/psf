INCLUDE_DIRECTORIES(
    ${PSF_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/testdata.h.in
        ${CMAKE_CURRENT_SOURCE_DIR}/testdata.h
        @ONLY IMMEDIATE
    )

#ADD_SUBDIRECTORY(testdata)

#### Sources
SET(SRCS_SPECTRUMALGORITHM SpectrumAlgorithm-test.cpp)
SET(SRCS_PEAKPARAMETER PeakParameter-test.cpp)
SET(SRCS_PEAKSHAPE PeakShape-test.cpp)
SET(SRCS_PEAKSHAPEFUNCTION  PeakShapeFunction-test.cpp)

MACRO(ADD_PSF_TEST name exe src)
    STRING(REGEX REPLACE "test_([^ ]+).*" "\\1" test "${exe}" )

    #build the test
    ADD_EXECUTABLE(${exe} ${src})
    #link the test
    TARGET_LINK_LIBRARIES(${exe} ms++ psf)
    
    #add test to global list of unit test
    ADD_TEST(${name} ${exe})

    #Add target for the test
    STRING(REGEX REPLACE "test_([^ ]+).*" "unit_\\1" unittest_target "${exe}" )
    ADD_CUSTOM_TARGET(${unittest_target} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${exe})
ENDMACRO(ADD_PSF_TEST name exe src)


#### Unit tests
ADD_PSF_TEST("PeakParameter" test_peakparameter ${SRCS_PEAKPARAMETER})
ADD_PSF_TEST("PeakShape" test_peakshape ${SRCS_PEAKSHAPE})
ADD_PSF_TEST("PeakShapeFunction" test_peakshapefunction ${SRCS_PEAKSHAPEFUNCTION})
ADD_PSF_TEST("SpectrumAlgorithm" test_spectrumalgorithm ${SRCS_SPECTRUMALGORITHM})

