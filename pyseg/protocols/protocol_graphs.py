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
import shutil
from glob import glob

from pwem.protocols import EMProtocol
from pyseg.convert import splitPysegStarFile
from pyseg.utils import createStarDirectories, genOutSplitStarFileName
from pyworkflow import BETA
from pyworkflow.protocol import FloatParam, PointerParam, LEVEL_ADVANCED, BooleanParam
from pyworkflow.utils import Message, makePath, moveFile
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
    def __init__(self):
        EMProtocol.__init__(self)
        ProtTomoBase.__init__(self)
        ProtTomoImportAcquisition.__init__(self)
        self._outStarDir = None
        self._inStarDir = None

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSegProt', PointerParam,
                      pointerClass='ProtPySegPreSegParticles',
                      label='Pre segmentation',
                      important=True,
                      allowsNull=False,
                      help='Pointer to preseg protocol.')
        form.addParam('keepOnlyReqFiles', BooleanParam,
                      label='Keep only required files?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to No, all the intermediate Disperse program resulting directories '
                           'will be kept in the extra folder.')

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
                       help='Vertex density within membranes. It allows to adjust simplification '
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
        starFileList = self._insertFunctionStep(self.convertInputStep)
        for starFile in starFileList:
            self._insertFunctionStep(self.pysegGraphs, starFile)
        self._insertFunctionStep(self.removeUnusedFilesStep)

    def convertInputStep(self):
        # Generate directories for input and output star files
        # Split the input file into n files, one per vesicle
        self._outStarDir, self._inStarDir = createStarDirectories(self._getExtraPath())
        return splitPysegStarFile(self._getPreSegStarFile(), self._inStarDir)

    def pysegGraphs(self, starFile):
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getGraphsCommand(starFile))
        # Fils returns the same star file name, so it will be renamed to avoid overwriting
        moveFile(self._getExtraPath('fil_mb_sources_to_no_mb_targets_net.star'),
                 genOutSplitStarFileName(self._outStarDir, starFile))

    def removeUnusedFilesStep(self):
        # Remove Disperse program intermediate result directories if requested
        if self.keepOnlyReqFiles.get():
            disperseDirs = glob(self._getExtraPath('disperse_*'))
            [shutil.rmtree(disperseDir) for disperseDir in disperseDirs]

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summaryMsg = []
        if self.isFinished():
            summaryMsg.append('Graphs were correctly generated.')

    # --------------------------- UTIL functions -----------------------------------
    def _getGraphsCommand(self, starFile):
        pixSize = self._getSamplingRate()/10
        graphsCmd = ' '
        graphsCmd += '%s ' % Plugin.getHome(GRAPHS_SCRIPT)
        graphsCmd += '--inStar %s ' % starFile
        graphsCmd += '--outDir %s ' % self._outStarDir
        graphsCmd += '--pixelSize %s ' % pixSize  # PySeg requires it in nm
        graphsCmd += '--sSig %s ' % self.sSig.get()
        graphsCmd += '--vDen %s ' % self.vDen.get()
        graphsCmd += '--veRatio %s ' % self.vRatio.get()
        graphsCmd += '--maxLen %s ' % self.maxLen.get()
        graphsCmd += '-j %s ' % self.numberOfThreads.get()
        return graphsCmd

    def _getPreSegStarFile(self):
        return self.inSegProt.get().getPresegOutputFile(self.inSegProt.get().getVesiclesCenteredStarFile())

    def _getSamplingRate(self):
        return self.inSegProt.get().outputSetofSubTomograms.getSamplingRate()
