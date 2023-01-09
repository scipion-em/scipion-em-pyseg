# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import glob
from os.path import abspath, join, basename
import numpy as np
from pwem.convert import transformations
from pwem.emlib.image import ImageHandler
from pyseg.constants import IN_STARS_DIR, OUT_STARS_DIR
from pyworkflow.utils import replaceExt, getExt, makePath, createLink
from reliontomo.constants import SHIFTX, SHIFTY, SHIFTZ, TILT, PSI, ROT

COMP_EXT_MASK_LIST = ['.mrc', '.em', '.rec']


def encodePresegArea(areaIndex):
    """EnumParam in form goes from 0 to 2 for the preseg areas codification, while
    Pyseg expects values higher or equal to 1, because 0 is used for the backgrounf"""
    return int(areaIndex) + 1


def getFinalMaskFileName(inMask):
    if getExt(inMask.getFileName()) in COMP_EXT_MASK_LIST:
        return abspath(inMask.getFileName())
    else:
        return abspath(replaceExt(inMask.getFileName(), '.mrc'))


def checkMaskFormat(inMask):
    if getExt(inMask.getFileName()) not in COMP_EXT_MASK_LIST:
        ih = ImageHandler()
        ih.convert(inMask.getFileName(), getFinalMaskFileName(inMask))


def createStarDirectories(extraPath):
    """Create a folder in extra directory to store the star files (one per vesicle) in which the input one will
    # be split and the same for the output star files"""
    splitStarDir = join(extraPath, IN_STARS_DIR)
    # The same for the output star files
    outStarDir = join(extraPath, OUT_STARS_DIR)
    makePath(splitStarDir, outStarDir)
    return outStarDir, splitStarDir


def getPrevPysegProtOutStarFiles(inDir, outDir):
    """inDir is the directory of the previous protocol output star files directory, while outDir is the input star
    files directory of the current protocol"""
    inStarList = glob.glob(join(inDir, '*.star'))
    outStarFiles = []
    for inStarFile in inStarList:
        outStarFile = join(outDir, basename(inStarFile))
        createLink(abspath(inStarFile), abspath(outStarFile))
        outStarFiles.append(outStarFile)

    return outStarFiles


def genOutSplitStarFileName(outDir, starFile):
    return join(outDir, basename(starFile))


# TODO: remove this once reliontomo3 is deprecated and import this method from reliontomo4 utils
def manageDims(fileName, z, n):
    if fileName.endswith('.mrc') or fileName.endswith('.map'):
        if z == 1 and n != 1:
            zDim = n
        else:
            zDim = z
    else:
        zDim = z

    return zDim




