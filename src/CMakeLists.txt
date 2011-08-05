SET(SRCS 
    BoxPeakShape.cpp
    ConstantModel.cpp
    GaussianPeakShape.cpp
    LinearSqrtModel.cpp
    LorentzianPeakShape.cpp
    PeakShapeFunction.cpp
    QuadraticModel.cpp
    SqrtModel.cpp
)

ADD_LIBRARY(psf ${SRCS})


