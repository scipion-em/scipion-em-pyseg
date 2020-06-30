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
    # Atfs = tfs.euler_matrix(np.deg2rad(rot), np.deg2rad(tilt), np.deg2rad(psi), 'szyz')
    A = eulerAngles2matrix(rot, tilt, psi, shiftx, shifty, shiftz)
    Ainv = np.linalg.inv(A)
    # A[0, 3] = shiftx
    # A[1, 3] = shifty
    # A[2, 3] = shiftz
    transform = Transform()
    transform.setMatrix(Ainv)
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


def eulerAngles2matrix(alpha, beta, gamma, shiftx, shifty, shiftz):
    A = np.empty([4, 4])
    A.fill(2)
    A[3, 3] = 1
    A[3, 0:3] = 0
    A[0, 3] = float(shiftx)
    A[1, 3] = float(shifty)
    A[2, 3] = float(shiftz)
    alpha = np.deg2rad(float(alpha))
    beta = np.deg2rad(float(beta))
    gamma = np.deg2rad(float(gamma))
    sa = np.sin(alpha)
    ca = np.cos(alpha)
    sb = np.sin(beta)
    cb = np.cos(beta)
    sg = np.sin(gamma)
    cg = np.cos(gamma)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    A[0, 0] = cg * cc - sg * sa
    A[0, 1] = cg * cs + sg * ca
    A[0, 2] = -cg * sb
    A[1, 0] = -sg * cc - cg * sa
    A[1, 1] = -sg * cs + cg * ca
    A[1, 2] = sg * sb
    A[2, 0] = sc
    A[2, 1] = ss
    A[2, 2] = cb

    return A