/*$Id$*/

/*
 * ModelMatrix.h
 *
 * Copyright (c) 2008 Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
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

