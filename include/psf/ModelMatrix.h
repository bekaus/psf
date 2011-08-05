#ifndef __MODELMATRIX_H__
#define __MODELMATRIX_H__
#include <ms++/config.h>

#include <vigra/matrix.hxx>

/**
 * Typedef to keep the ms++ interface clean of VIGRA code
 */
namespace psf
{
typedef vigra::linalg::Matrix<double> ModelMatrix;
}

#endif /*__MODELMATRIX_H__*/

