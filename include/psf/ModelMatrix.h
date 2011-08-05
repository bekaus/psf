#ifndef __MODELMATRIX_H__
#define __MODELMATRIX_H__
#include <psf/config.h>

#include <vigra/matrix.hxx>

/**
 * Typedef to keep the psf interface clean of VIGRA code
 */
namespace psf
{
typedef vigra::linalg::Matrix<double> ModelMatrix;
}

#endif /*__MODELMATRIX_H__*/

