SET(SRCS-MSPP 
    SparseSpectrum.cpp
    PeakShapeFunction.cpp
    OrbiBoxPeakShapeFunction.cpp
    OrbiBoxPeakShapeFunction.cpp
)

SET(SRCS 
    BoxPeakShape.cpp
    ConstantModel.cpp
    GaussianPeakShape.cpp
    LinearSqrtModel.cpp
    LorentzianPeakShape.cpp
    QuadraticModel.cpp
    SpectralPeak.cpp
    SpectrumAlgorithm.cpp
    SqrtModel.cpp
)

ADD_LIBRARY(ms++ ${SRCS-MSPP})
ADD_LIBRARY(psf ${SRCS})

TARGET_LINK_LIBRARIES(psf ms++)


