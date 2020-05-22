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
import math
import numpy as np
from pwem import Domain
from pwem.emlib.image.image_handler import ImageHandler
from pwem.objects.data import Transform
import pwem.convert.transformations as tfs
Coordinate3D = Domain.importFromPlugin("tomo.objects", "Coordinate3D")
TomoAcquisition = Domain.importFromPlugin("tomo.objects", "TomoAcquisition")

def readStarfileRow(self, item):
    nline = next(self.fhTable)
    nline = nline.rstrip()

    item.setFileName(nline.split()[2])

    coordinate3d = Coordinate3D()
    volName = int(nline.split()[0])
    ctf3d = nline.split()[1]
    x = nline.split()[3]
    y = nline.split()[4]
    z = nline.split()[5]
    coordinate3d.setVolName(volName)
    coordinate3d.setX(float(x))
    coordinate3d.setY(float(y))
    coordinate3d.setZ(float(z))
    # Extended attribute in prot import starfile
    coordinate3d.setCtf3D(ctf3d)
    item.setCoordinate3D(coordinate3d)

    shiftx = nline.split()[12]
    shifty = nline.split()[13]
    shiftz = nline.split()[14]
    tilt = nline.split()[6]
    psi = nline.split()[8]
    rot = nline.split()[10]
    # Check angles order in euler_matrix
    A = tfs.euler_matrix(tilt, psi, rot, shiftx, shifty, shiftz)
    transform = Transform()



