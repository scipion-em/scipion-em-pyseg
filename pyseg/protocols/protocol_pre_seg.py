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

from os.path import abspath

from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FileParam, NumericListParam, IntParam, FloatParam, GT, LEVEL_ADVANCED, EnumParam, \
    PointerParam
from pyworkflow.utils import Message, removeBaseExt, removeExt
from scipion.constants import PYTHON

from pyseg import Plugin
from pyseg.constants import PRESEG_SCRIPT, TOMOGRAM, PYSEG_LABEL, VESICLE, NOT_FOUND, \
    PYSEG_OFFSET_X, PYSEG_OFFSET_Y, PYSEG_OFFSET_Z, MASK, RLN_ORIGIN_X, RLN_ORIGIN_Y, \
    RLN_ORIGIN_Z
from relion.convert import Table

SEG_FROM_SCIPION = 0
SEG_FROM_STAR = 1


class ProtPySegPreSegParticles(EMProtocol):
    """pre-process segmented circular membranes"""

    _label = 'preseg circular membranes'
    _devStatus = BETA
    _starFile = None

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('segmentationFrom', EnumParam,
                      choices=['Scipion Protocol', 'Star file'],
                      default=0,
                      label='Choose pre seg data source',
                      important=True,
                      display=EnumParam.DISPLAY_HLIST)
        form.addParam('inTomoMasks', PointerParam,
                      pointerClass='SetOfTomoMasks',
                      label='Segmented and annotated tomograms',
                      condition='segmentationFrom == %i' % SEG_FROM_SCIPION,
                      important=True,
                      allowsNull=False,
                      help='Pointer to segmented and annotated tomograms within Scipion.')
        form.addParam('inStar', FileParam,
                      label='Seg particles star file',
                      condition='segmentationFrom == %i' % SEG_FROM_STAR,
                      important=True,
                      allowsNull=False,
                      help='Star file obtained in PySeg graphs step.')
        group = form.addGroup('Sub-volume splitting')
        group.addParam('spSplit', NumericListParam,
                       default='-1',
                       allowsNull=False,
                       label='Number of splits (X, Y, Z)',
                       help='Parts in which the tomogram will be split, respecting X, Y and Z axis. Value -1'
                            'is used to indicate no splitting.')
        group.addParam('spOffVoxels', IntParam,
                       label='Offset voxels',
                       allowsNull=False,
                       validators=[GT(0)],
                       help='Margin to ensure that the desired entities, e. g. membranes, proteins, are included.')
        group = form.addGroup('Membrane segmentation')
        group.addParam('sgVoxelSize', FloatParam,
                       label='Voxel size (Å/voxel)',
                       validators=[GT(0)],
                       important=True,
                       allowsNull=False)
        group.addParam('sgThreshold', IntParam,
                       default=-1,
                       label='Density threshold',
                       expertLevel=LEVEL_ADVANCED,
                       help='All the voxels with density equal to or higher than the threshold are set to 1. '
                            'The remaining voxels are set to 0.')
        group.addParam('sgSizeThreshold', IntParam,
                       default=-1,
                       label='Size threshold (voxels)',
                       condition='sgThreshold > 0',
                       help='It sets the minimal size for a component to be considered as membrane.')
        group.addParam('sgMembThk', FloatParam,
                       label='Segmented membrane thickness (Å)',
                       allowsNull=False,
                       validators=[GT(0)])
        group.addParam('sgMembNeigh', FloatParam,
                       label='Segmented mebmrane neighbours (Å)',
                       allowsNull=False,
                       validators=[GT(0)],
                       help='Thickness around the membrane to represent the in-membrane and out-membrane surroundings '
                            'desired to be included in the analysis.')

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('pysegPreSegStep')
        self._insertFunctionStep('getMembraneCenterStep')
        self._insertFunctionStep('pysegPreSegCenteredStep')

    def _initialize(self):
        if self.segmentationFrom.get() == SEG_FROM_SCIPION:
            self._starFile = self._convertInput2Star()
        else:
            self._starFile = self.inStar.get()

    def pysegPreSegStep(self):
        outDir = self._getExtraPath()
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPreSegCmd(self._starFile, outDir))

    def getMembraneCenterStep(self):
        inStar = self.getPresegOutputFile(self._starFile)
        self._findVesicleCenter(self._starFile, inStar)

    def pysegPreSegCenteredStep(self):
        inStar = abspath(self.getVesiclesCenteredStarFile())
        outDir = self._getExtraPath()

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPreSegCmd(inStar, outDir))

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        # if self.isFinished():
        #     summary.append('Generated files location:\n'
        #                    'Subtomograms files directory: %s\n'
        #                    'Star file: %s' %
        #                    (self._getExtraPath(POST_REC_OUT),
        #                     self._getExtraPath(POST_REC_OUT + '.star')))
        return summary

    # --------------------------- UTIL functions -----------------------------------

    def _getPreSegCmd(self, inStar, outDir):
        preSegCmd = ' '
        preSegCmd += '%s ' % Plugin.getHome(PRESEG_SCRIPT)
        preSegCmd += '--inStar %s ' % inStar
        preSegCmd += '--outDir %s ' % outDir
        preSegCmd += '--spSplit %s ' % self.spSplit.get()
        preSegCmd += '--spOffVoxels %s ' % self.spOffVoxels.get()
        preSegCmd += '--sgVoxelSize %s ' % (int(self.sgVoxelSize.get())/10)  # required in nm
        preSegCmd += '--sgThreshold %s ' % self.sgThreshold.get()
        preSegCmd += '--sgSizeThreshold %s ' % self.sgSizeThreshold.get()
        preSegCmd += '--sgMembThk %s ' % self._checkValue4PySeg(self.sgMembThk.get())  # required in nm
        preSegCmd += '--sgMembNeigh %s ' % self._checkValue4PySeg(self.sgMembNeigh.get())  # required in nm

        return preSegCmd

    def _convertInput2Star(self):
        inStar = self._getExtraPath('inStar.star')
        outputTable = self._createTable(isConvertingInput=True)

        # Generate the star file with the vesicles centered for the second pre_seg execution
        for tomoMask in self.inTomoMasks.get():
            vesicle = tomoMask.getFileName()

            for materialIndex in self._getMaterialsList(vesicle):   # Get annotated materials from txt file and add one
                # line for each one

                # Add row to output table
                outputTable.addRow(tomoMask.getVolName(),
                                   vesicle,
                                   int(materialIndex),
                                   vesicle
                                   )

        outputTable.write(inStar)
        return inStar

    @staticmethod
    def _getMaterialsList(vesicle):
        # Get annotated materials from txt file and add one line for each one
        materialsFile = removeExt(vesicle) + '.txt'
        with open(materialsFile) as matFile:
            materialsList = matFile.read()

        # Expected format is a string like 'ind1,ind2,...,indn\n, so it's necessary to transform it into a list of
        # material indices
        return materialsList.replace('\n', '').split(',')

    def _findVesicleCenter(self, starFileInit, starFilePreseg1):
        ih = ImageHandler()
        outputTable = self._createTable()
        # Read preseg (vesicles not centered) table
        presegTable = Table()
        presegTable.read(starFilePreseg1)
        # Read initial data
        initTable = Table()
        initTable.read(starFileInit)

        # Generate the star file with the vesicles centered for the second pre_seg execution
        for row, rowp in zip(initTable, presegTable):
            tomo = row.get(TOMOGRAM, NOT_FOUND)
            vesicle = row.get(VESICLE, NOT_FOUND)
            materialIndex = row.get(PYSEG_LABEL, NOT_FOUND)

            # Get the upper left corner of the vesicle in the original tomogram: _psSegOffX #7 _psSegOffY #8 _psSegOffZ
            # #9 from output star file of pre_tomos_seg.py that do not distinguish inner and outer densities
            xdimCorner = rowp.get(PYSEG_OFFSET_X, 0)
            ydimCorner = rowp.get(PYSEG_OFFSET_Y, 0)
            zdimCorner = rowp.get(PYSEG_OFFSET_Z, 0)

            # Get the box dimensions, be sure that the Dimensions format is Dimensions: 239 x 1 x 298 x 298
            # ((N)Objects x (Z)Slices x (Y)Rows x (X)Columns)
            x, y, z, _ = ih.getDimensions(rowp.get(MASK, NOT_FOUND))

            # Add row to output table
            outputTable.addRow(tomo,
                               vesicle,
                               materialIndex,
                               vesicle,
                               xdimCorner + x / 2,
                               ydimCorner + y / 2,
                               zdimCorner + z / 2,
                               )

        outputTable.write(self.getVesiclesCenteredStarFile())

    def getPresegOutputFile(self, inStar):
        return self._getExtraPath(removeBaseExt(inStar) + '_pre.star')

    def getVesiclesCenteredStarFile(self):
        return self._getExtraPath('presegVesiclesCentered.star')

    @ staticmethod
    def _createTable(isConvertingInput=False):
        if isConvertingInput:
            # Headers for the input star file to the first pre seg when the input data is from Scipion
            return Table(columns=[TOMOGRAM,
                                  VESICLE,
                                  PYSEG_LABEL,
                                  MASK])
        else:
            # Headers for pySeg pre_seg centered star file
            return Table(columns=[TOMOGRAM,
                                  VESICLE,
                                  PYSEG_LABEL,
                                  MASK,
                                  RLN_ORIGIN_X,
                                  RLN_ORIGIN_Y,
                                  RLN_ORIGIN_Z
                                  ])

    @staticmethod
    def _checkValue4PySeg(value):
        return value if value == -1 else value/10
