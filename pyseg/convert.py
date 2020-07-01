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
    A = eulerAngles2matrix(rot, tilt, psi)
    As = [float(shiftx), float(shifty), float(shiftz)]
    A = np.column_stack((A, As))
    A0 = [0, 0, 0, 1]
    A = np.vstack((A, A0))
    transform = Transform()
    transform.setMatrix(A)
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


def eulerAngles2matrix(alpha, beta, gamma):
    A = np.empty([3, 3])
    A.fill(0)

    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)

    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sg = np.sin(gamma)
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


def matrix2eulerAngles(A):
    epsilon = 16*np.exp(-5)

    if np.abs(A[1, 1]) > epsilon:
        abs_sb = np.sqrt((-A[2, 2] * A[1, 2] * A[2, 1] - A[0, 2] * A[2, 0]) / A[1, 1])

    elif np.abs(A[0, 1]) > epsilon:
        abs_sb = np.sqrt((-A[2, 1] * A[2, 2] * A[0, 2] + A[2, 0] * A[1, 2]) / A[0, 1])
    elif np.abs(A[0, 0]) > epsilon:
        abs_sb = np.sqrt((-A[2, 0] * A[2, 2] * A[0, 2] - A[2, 1] * A[1, 2]) / A[0, 0])
    else:
        exit(1)

    if abs_sb > epsilon:
        beta = np.arctan2(abs_sb, A[2, 2])
        alpha = np.arctan2(A[2, 1] / abs_sb, A[2, 0] / abs_sb)
        gamma = np.arctan2(A[1, 2] / abs_sb, -A[0, 2] / abs_sb)

    else:
        alpha = 0
        beta = 0
        gamma = np.arctan2(A[1, 0], A[0, 0])

    gamma = np.rad2deg(gamma)
    beta = np.rad2deg(beta)
    alpha = np.rad2deg(alpha)

    return gamma, beta, alpha
