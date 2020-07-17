# **************************************************************************
# *
# * Authors:  Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import numpy as np
from pwem.objects.data import Transform, String
import pwem.convert.transformations as tfs
from relion.convert import Table
from tomo.objects import SubTomogram, Coordinate3D, TomoAcquisition

RELION_TOMO_LABELS = ['rlnMicrographName',
                          'rlnCoordinateX',
                          'rlnCoordinateY',
                          'rlnCoordinateZ',
                          'rlnImageName',
                          'rlnCtfImage',
                          'rlnMagnification',
                          'rlnDetectorPixelSize',
                          'rlnAngleRot',
                          'rlnAngleTilt',
                          'rlnAnglePsi',
                          'rlnOriginX',
                          'rlnOriginY',
                          'rlnOriginZ',
                          ]
FILE_NOT_FOUND = 'file_not_found'


def readStarFile(starFile, outputSubTomogramsSet, invert=True):
    warningMsg = ''
    tomoTable = Table()
    tomoTable.read(starFile)
    if not tomoTable.hasAllColumns(RELION_TOMO_LABELS):
        missingCols = [name for name in RELION_TOMO_LABELS if name not in tomoTable.getColumnNames()]
        warningMsg = 'Columns %s\nwere not found in the star file provided.\nThe corresponding numerical ' \
                     'values will be considered as 0.' \
                     % '  '.join(['*' + colName + '*' for colName in missingCols])

    for counter, row in enumerate(tomoTable):
        subtomo = SubTomogram()
        coordinate3d = Coordinate3D()
        transform = Transform()

        volname = row.get('rlnMicrographName', FILE_NOT_FOUND)
        filename = row.rlnImageName
        x = row.get('rlnCoordinateX', 0)
        y = row.get('rlnCoordinateY', 0)
        z = row.get('rlnCoordinateZ', 0)
        ctf3d = row.get('rlnCtfImage', FILE_NOT_FOUND)
        coordinate3d.setX(float(x))
        coordinate3d.setY(float(y))
        coordinate3d.setZ(float(z))
        coordinate3d._3dcftMrcFile = String(ctf3d) # This extended attribute is used for the ctf3d generation in Relion
        shiftx = row.get('rlnOriginX', 0)
        shifty = row.get('rlnOriginY', 0)
        shiftz = row.get('rlnOriginZ', 0)
        tilt = row.get('rlnAngleTilt', 0)
        psi = row.get('rlnAnglePsi', 0)
        rot = row.get('rlnAngleRot', 0)
        shifts = (float(shiftx), float(shifty), float(shiftz))
        angles = (float(rot), float(tilt), float(psi))
        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        if invert:
            M[0, 3] = -shifts[0]
            M[1, 3] = -shifts[1]
            M[2, 3] = -shifts[2]
            M = np.linalg.inv(M)
        else:
            M[0, 3] = shifts[0]
            M[1, 3] = shifts[1]
            M[2, 3] = shifts[2]
        transform.setMatrix(M)

        subtomo.setVolName(volname)
        subtomo.setFileName(filename)
        subtomo.setCoordinate3D(coordinate3d)
        subtomo.setTransform(transform)
        subtomo.setAcquisition(TomoAcquisition())
        subtomo.setClassId(row.get('rlnClassNumber', 0))
        outputSubTomogramsSet.append(subtomo)

    return warningMsg

