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

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FloatParam, EnumParam, PointerParam, FileParam, LEVEL_ADVANCED
from pyworkflow.utils import Message, makePath
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition

from pyseg import Plugin
from pyseg.constants import GRAPHS_SCRIPT, FROM_SCIPION, FROM_STAR_FILE


class ProtPySegGraphs(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """analyze a GraphMCF (Mean Cumulative Function) from a segmented membrane"""

    _label = 'graphs'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('presegFrom', EnumParam,
                      choices=['Scipion Protocol', 'Star file'],
                      default=FROM_SCIPION,
                      label='Choose preSeg data source',
                      important=True,
                      display=EnumParam.DISPLAY_HLIST)
        form.addParam('inSegProt', PointerParam,
                      pointerClass='ProtPySegPreSegParticles',
                      label='Pre segmentation',
                      condition='presegFrom == %s' % FROM_SCIPION,
                      important=True,
                      allowsNull=False,
                      help='Pointer to preseg protocol.')
        form.addParam('inStar', FileParam,
                      label='Seg particles star file',
                      condition='presegFrom == %s' % FROM_STAR_FILE,
                      important=True,
                      allowsNull=False,
                      help='Star file obtained in PySeg seg step.')
        form.addParam('pixelSize', FloatParam,
                      label='Pixel size (Å/voxel)',
                      important=True,
                      condition='presegFrom == %s' % FROM_STAR_FILE,
                      help='Input tomograms voxel size (Å/voxel)')

        group = form.addGroup('Graphs parameters', expertLevel=LEVEL_ADVANCED)
        group.addParam('sSig', FloatParam,
                       label='Sigma for gaussian filtering',
                       default=1,
                       allowsNull=False,
                       expertLevel=LEVEL_ADVANCED,
                       help='Sigma for Gaussian fltering input tomograms. It allows to smooth '
                            'small and irrelevant features and increases SNR.')
        group.addParam('vDen', FloatParam,
                       label='Vertex density within membranes (nm³)',
                       default=0.0035,
                       allowsNull=False,
                       expertLevel=LEVEL_ADVANCED,
                       help='Vertex density within membranes. It allows to adjust simplifcation '
                            'adaptively for every tomogram.')
        group.addParam('vRatio', FloatParam,
                       label='Avg ratio vertex/edge of graph within membrane',
                       default=4,
                       allowsNull=False,
                       expertLevel=LEVEL_ADVANCED,
                       help='Averaged ratio vertex/edge in the graph within membrane.')
        group.addParam('maxLen', FloatParam,
                       label='Shortest distance to membrane (nm)',
                       default=10,
                       allowsNull=False,
                       expertLevel=LEVEL_ADVANCED,
                       help='Maximum euclidean shortest distance to membrane in nm.')

        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.pysegGraphs.__name__)

    def pysegGraphs(self):
        # Generate output dir
        outDir = self._getExtraPath()
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getGraphsCommand(outDir))

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summaryMsg = []
        if self.isFinished():
            summaryMsg.append('Graphs were correctly generated.')

    def _validate(self):
        validationMsg = []
        if self.presegFrom.get() == FROM_STAR_FILE:
            voxelSize = self.pixelSize.get()
            msg = 'Pixel size must be greater than 0.'
            if voxelSize:
                if voxelSize <= 0:
                    validationMsg.append(msg)
            else:
                validationMsg.append(msg)

        return validationMsg

    # --------------------------- UTIL functions -----------------------------------
    def _getGraphsCommand(self, outDir):
        pixSize = self._getSamplingRate()/10
        graphsCmd = ' '
        graphsCmd += '%s ' % Plugin.getHome(GRAPHS_SCRIPT)
        graphsCmd += '--inStar %s ' % self._getPreSegStarFile()  # self.inStar.get()
        graphsCmd += '--outDir %s ' % outDir
        graphsCmd += '--pixelSize %s ' % pixSize  # PySeg requires it in nm
        graphsCmd += '--sSig %s ' % self.sSig.get()
        graphsCmd += '--vDen %s ' % self.vDen.get()
        graphsCmd += '--veRatio %s ' % self.vRatio.get()
        graphsCmd += '--maxLen %s ' % self.maxLen.get()
        graphsCmd += '-j %s ' % self.numberOfThreads.get()
        return graphsCmd

    def _getPreSegStarFile(self):
        if self.presegFrom.get() == FROM_SCIPION:
            return self.inSegProt.get().getPresegOutputFile(self.inSegProt.get().getVesiclesCenteredStarFile())
        else:
            return self.inStar.get()

    def _getSamplingRate(self):
        if self.presegFrom.get() == FROM_SCIPION:
            return self.inSegProt.get().outputSetofSubTomograms.getSamplingRate()
        else:
            return self.pixelSize.get()
