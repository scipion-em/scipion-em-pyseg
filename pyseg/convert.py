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
from os.path import dirname, abspath, isabs

import numpy as np
from pwem.emlib.image import ImageHandler
from pwem.objects.data import Transform, String, Integer
import pwem.convert.transformations as tfs
from os.path import join
from pyworkflow.utils import createAbsLink
from relion.convert import Table
from tomo.objects import SubTomogram, Coordinate3D, TomoAcquisition


FILE_NOT_FOUND = 'file_not_found'

# Star file fields
TOMO_NAME = 'rlnMicrographName'
COORD_X = 'rlnCoordinateX'
COORD_Y = 'rlnCoordinateY'
COORD_Z = 'rlnCoordinateZ'
SUBTOMO_NAME = 'rlnImageName'
CTF_MISSING_WEDGE = 'rlnCtfImage'
MAGNIFICATION = 'rlnMagnification'
PIXEL_SIZE = 'rlnDetectorPixelSize'
ROT = 'rlnAngleRot'
TILT = 'rlnAngleTilt'
PSI = 'rlnAnglePsi'
SHIFTX = 'rlnOriginX'
SHIFTY = 'rlnOriginY'
SHIFTZ = 'rlnOriginZ'
PS_SEG_IMAGE = '_psSegImage'

RELION_TOMO_LABELS = [TOMO_NAME,
                      COORD_X,
                      COORD_Y,
                      COORD_Z,
                      SUBTOMO_NAME,
                      CTF_MISSING_WEDGE,
                      MAGNIFICATION,
                      PIXEL_SIZE,
                      ROT,
                      TILT,
                      PSI,
                      SHIFTX,
                      SHIFTY,
                      SHIFTZ]

PYSEG_PICKING_LABELS = [TOMO_NAME,
                        SUBTOMO_NAME,
                        PS_SEG_IMAGE,
                        COORD_X,
                        COORD_Y,
                        COORD_Z,
                        ROT,
                        TILT,
                        PSI]

# Star files coding
RELION_SUBTOMO_STAR = 0
PYSEG_PICKING_STAR = 1


def readStarFile(prot, outputSubTomogramsSet, fileType, starFile=None, invert=True):
    warningMsg = None
    # Star file can be provided by the user or not, depending on the protocol invoking this method
    if not starFile:
        starFile = prot.starFile.get()

    # If the star file is currently in the extra folder of another protocol execution, the paths
    # generated with method _getExtraPath() will become wrong, but it doesn't have to be located there
    if 'extra' in starFile:
        starPath = ''
        isExtra = True
    else:
        starPath = dirname(starFile) + '/'
        isExtra = False

    tomoTable = Table()
    tomoTable.read(starFile)

    if fileType == RELION_SUBTOMO_STAR:
        labels = RELION_TOMO_LABELS
        _processRelionTomoStar(prot, outputSubTomogramsSet, tomoTable, starPath, isExtra, invert)
    elif fileType == PYSEG_PICKING_STAR:
        labels = PYSEG_PICKING_LABELS

    if not tomoTable.hasAllColumns(labels):
        missingCols = [name for name in labels if name not in tomoTable.getColumnNames()]
        warningMsg = 'Columns %s\nwere not found in the star file provided.\nThe corresponding numerical ' \
                     'values will be considered as 0.' \
                     % '  '.join(['*' + colName + '*' for colName in missingCols])

    return warningMsg


def manageIhDims(fileName, z, n):
    if fileName.endswith('.mrc') or fileName.endswith('.map'):
        fileName += ':mrc'
        if z == 1 and n != 1:
            zDim = n
        else:
            zDim = z
    else:
        zDim = z

    return zDim, fileName


def genAbsLink(fileName, newFileName):
    if fileName.endswith(':mrc'):
        fileName = fileName[:-4]
    createAbsLink(fileName, newFileName)


def _processRelionTomoStar(prot, outputSubTomogramsSet, tomoTable, starPath, isExtra, invert):
    ih = ImageHandler()
    samplingRate = outputSubTomogramsSet.getSamplingRate()
    for counter, row in enumerate(tomoTable):
        subtomo = SubTomogram()
        coordinate3d = Coordinate3D()
        transform = Transform()
        origin = Transform()

        volname = join(starPath, row.get(TOMO_NAME, FILE_NOT_FOUND))
        subtomoFn = row.get(SUBTOMO_NAME, FILE_NOT_FOUND)
        subtomoAbsFn = join(starPath, subtomoFn)
        if not isExtra:
            subtomoAbsFn = prot._getExtraPath(subtomoAbsFn)
        if not isabs(subtomoAbsFn):
            subtomoAbsFn = abspath(subtomoAbsFn)
        x = row.get(COORD_X, 0)
        y = row.get(COORD_Y, 0)
        z = row.get(COORD_Z, 0)
        ctf3d = row.get(CTF_MISSING_WEDGE, FILE_NOT_FOUND)
        coordinate3d.setX(float(x))
        coordinate3d.setY(float(y))
        coordinate3d.setZ(float(z))
        coordinate3d._3dcftMrcFile = String(join(starPath, ctf3d))  # Used for the ctf3d generation in Relion
        shiftx = row.get(SHIFTX, 0)
        shifty = row.get(SHIFTY, 0)
        shiftz = row.get(SHIFTZ, 0)
        tilt = row.get(TILT, 0)
        psi = row.get(PSI, 0)
        rot = row.get(ROT, 0)
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
        subtomo.setCoordinate3D(coordinate3d)
        subtomo.setTransform(transform)
        subtomo.setAcquisition(TomoAcquisition())
        subtomo.setClassId(row.get('rlnClassNumber', 0))
        subtomo.setSamplingRate(samplingRate)

        # Make link if necessary (only when the star file is out of a Scipion results dir)
        uniqueSubtomoFn = subtomoAbsFn
        if not isExtra:
            uniqueSubtomoFn = subtomoFn.replace("/", "_").replace("..", "")
            genAbsLink(subtomoAbsFn, uniqueSubtomoFn)

        # Set the origin and the dimensions of the current subtomogram
        x, y, z, n = ih.getDimensions(subtomoAbsFn)
        zDim, filename = manageIhDims(uniqueSubtomoFn, z, n)
        origin.setShifts(x / -2. * samplingRate, y / -2. * samplingRate, zDim / -2. * samplingRate)
        subtomo.setOrigin(origin)

        subtomo.setFileName(filename)
        # if subtomo is in a vesicle
        if 'tid_' in subtomoFn:
            vesicleId = subtomoFn.split('tid_')[1]
            vesicleId = vesicleId[0]
            scoor = subtomo.getCoordinate3D()
            scoor.setGroupId(vesicleId)
            subtomo.setCoordinate3D(scoor)

        # Add current subtomogram to the output set
        outputSubTomogramsSet.append(subtomo)

