
INCLUDE_DIRECTORIES(
    ${MSPP_SOURCE_DIR}/src/parsers/mgf
)

SET(SRCS 
    BoxPeakShape.cpp
    ConstantModel.cpp
    GaussianPeakShape.cpp
    LinearSqrtModel.cpp
    LorentzianPeakShape.cpp
    OrbiBoxPeakShapeFunction.cpp
    PeakShapeFunction.cpp
    QuadraticModel.cpp
    SparseSpectrum.cpp
    SpectralPeak.cpp
    SpectrumAlgorithm.cpp
    SqrtModel.cpp
)

#for MSVC to define a "Header Files" subfolder:
SOURCE_GROUP("Header Files" FILES mspp_header)

ADD_LIBRARY(ms++ ${SRCS} ${mspp_header})

#TARGET_LINK_LIBRARIES(ms++)

#TARGET_LINK_LIBRARIES(ms++ ${IMPEX_LIB} ${PNG_LIBRARY})

FILE(GLOB headers "${MSPP_SOURCE_DIR}/include/ms++/*.h")
INSTALL(FILES ${headers} DESTINATION include/ms++)

#install the library in lib/ under unix and in bin/ under windows
INSTALL(TARGETS ms++
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
)
