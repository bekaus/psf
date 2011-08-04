SET(SRCS-MSPP 
    SparseSpectrum.cpp
)

SET(SRCS 
    BoxPeakShape.cpp
    ConstantModel.cpp
    GaussianPeakShape.cpp
    LinearSqrtModel.cpp
    LorentzianPeakShape.cpp
    PeakShapeFunction.cpp
    QuadraticModel.cpp
    SpectrumAlgorithm.cpp
    SqrtModel.cpp
)

ADD_LIBRARY(ms++ ${SRCS-MSPP})
ADD_LIBRARY(psf ${SRCS})

TARGET_LINK_LIBRARIES(psf ms++)


