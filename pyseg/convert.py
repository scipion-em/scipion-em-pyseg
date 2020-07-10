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
from os.path import join
import numpy as np
from pwem import Domain
from pwem.objects.data import Transform, String
import pwem.convert.transformations as tfs

Coordinate3D = Domain.importFromPlugin("tomo.objects", "Coordinate3D")
TomoAcquisition = Domain.importFromPlugin("tomo.objects", "TomoAcquisition")


def readStarfileRow(nline, item, path, headerDict):
    nline = nline.rstrip()
    volname = join(path, nline.split()[headerDict.get('_rlnMicrographName')])
    item.setVolName(volname)
    filename = join(path, nline.split()[headerDict.get('_rlnImageName')])
    item.setFileName(filename)
    coordinate3d = Coordinate3D()
    x = nline.split()[headerDict.get('_rlnCoordinateX')]
    y = nline.split()[headerDict.get('_rlnCoordinateY')]
    z = nline.split()[headerDict.get('_rlnCoordinateZ')]
    coordinate3d.setX(float(x))
    coordinate3d.setY(float(y))
    coordinate3d.setZ(float(z))
    ctf3d = nline.split()[headerDict.get('_rlnCtfImage')]
    # This extended attribute should match with ctf3d generation in Relion
    coordinate3d._3dcftMrcFile = String(join(path, ctf3d))
    item.setCoordinate3D(coordinate3d)
    shiftx = float(nline.split()[headerDict.get('_rlnOriginX')])
    shifty = float(nline.split()[headerDict.get('_rlnOriginY')])
    shiftz = float(nline.split()[headerDict.get('_rlnOriginZ')])
    tilt = float(nline.split()[headerDict.get('_rlnAngleTilt')])
    psi = float(nline.split()[headerDict.get('_rlnAnglePsi')])
    rot = float(nline.split()[headerDict.get('_rlnAngleRot')])
    shifts = (shiftx, shifty, shiftz)
    angles = (rot, tilt, psi)
    invert = True
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
    transform = Transform()
    transform.setMatrix(M)
    item.setTransform(transform)
    item.setClassId(int(nline.split()[headerDict.get('_rlnClassNumber')]))
    acq = TomoAcquisition()
    item.setAcquisition(acq)


def readStarfileHeader(fhStar):
    w1 = True
    while(w1):  # Read lines until they contain first column name
        line = next(fhStar)
        if line.startswith('loop_'):
            w1 = False
            break

    w2 = True
    headerList = []
    while(w2):  # Read all lines with column names and store the names in a list
        line = next(fhStar)
        if line.startswith('_rln'):
            headerList.append(line)
        else:
            w2 = False
            break

    headerDict = {}
    for i, colName in enumerate(headerList):
        headerDict[colName.split()[0]] = i

    return headerDict, line
