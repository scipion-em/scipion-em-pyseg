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
from numpy import deg2rad
from pwem import Domain
from pwem.objects.data import Transform, String
import pwem.convert.transformations as tfs
Coordinate3D = Domain.importFromPlugin("tomo.objects", "Coordinate3D")
TomoAcquisition = Domain.importFromPlugin("tomo.objects", "TomoAcquisition")


def readStarfileRow(nline, item, path):
    nline = nline.rstrip()
    volname = join(path, nline.split()[0])
    item.setVolName(volname)
    filename = join(path, nline.split()[2])
    item.setFileName(filename)
    coordinate3d = Coordinate3D()
    x = nline.split()[3]
    y = nline.split()[4]
    z = nline.split()[5]
    coordinate3d.setX(float(x))
    coordinate3d.setY(float(y))
    coordinate3d.setZ(float(z))
    ctf3d = nline.split()[1]
    # This extended attribute should match with ctf3d generation in Relion
    coordinate3d._3dcftMrcFile = String(join(path, ctf3d))
    item.setCoordinate3D(coordinate3d)
    shiftx = float(nline.split()[12])
    shifty = float(nline.split()[13])
    shiftz = float(nline.split()[14])
    tilt = float(nline.split()[6])
    psi = float(nline.split()[8])
    rot = float(nline.split()[10])
    A = tfs.euler_matrix(deg2rad(rot), deg2rad(tilt), deg2rad(psi), 'szyz')
    A[0, 3] = shiftx
    A[1, 3] = shifty
    A[2, 3] = shiftz
    transform = Transform()
    item.setTransform(transform)
    item.setClassId(int(nline.split()[15]))
    acq = TomoAcquisition()
    item.setAcquisition(acq)
