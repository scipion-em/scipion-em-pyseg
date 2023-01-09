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
from enum import Enum

from pwem.protocols import EMProtocol, PointerParam
from pyseg.convert import readPysegSubtomograms
from pyseg.utils import getFinalMaskFileName, checkMaskFormat
from pyworkflow import BETA
from pyworkflow.protocol import String, FloatParam, LE, GE
from pyworkflow.utils import Message, makePath
from reliontomo.convert import createWriterTomo
from scipion.constants import PYTHON
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_OUT, POST_REC_SCRIPT_MEMB_ATT, SEE_METHODS_TAB


class outputObjects(Enum):
    subtomograms = SetOfSubTomograms


class ProtPySegPostRecParticles(EMProtocol, ProtTomoBase):
    """post-process already reconstructed particles: rot angle randomization and membrane suppression"""

    _label = 'posrec'
    _devStatus = BETA
    inStarName = 'input_particles.star'
    warningMsg = None
    subtomoSet = None

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label="Input subtomograms",
                      help='Select the input subtomograms from the project.')
        form.addParam('inMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask',
                      important=True,
                      allowsNull=False,
                      help='Mask used for the post processing')
        group = form.addGroup('Membrane suppression')
        group.addParam('mbMask', PointerParam,
                       pointerClass='VolumeMask',
                       allowsNull=True,
                       label='Mask for membrane supression')
        group.addParam('mbSupFactor', FloatParam,
                       label='Membrane suppression factor [0, 1]',
                       default=0.5,
                       allowsNull=False,
                       validators=[GE(0), LE(1)],
                       condition='mbMask',
                       help='Value 0 suppress the area corresponding to the suppression mask, while higher values '
                            'up to 1 attenuate it.')
        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        outStar = self._getExtraPath(POST_REC_OUT + '.star')
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.pysegPostRec, outStar)
        self._insertFunctionStep(self.createOutputStep, outStar)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        # Check masks format and convert if necessary
        checkMaskFormat(self.inMask.get())  # Subtomogram mask
        if self.mbMask.get():
            checkMaskFormat(self.mbMask.get())  # Membrane mask for attenuation (optional)
        # Write star from set of subtomograms
        subtomoSet = self.inputSubtomos.get()
        subTomoStar = self._getExtraPath(self.inStarName)
        writer = createWriterTomo(isPyseg=True)
        writer.subtomograms2Star(subtomoSet, subTomoStar)
        # Convert the mask format if necessary
        checkMaskFormat(self.inMask.get())

    def pysegPostRec(self, outStar):
        # Generate output subtomo dir
        outDir = self._getExtraPath(POST_REC_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getCommand(outDir, outStar))

    def createOutputStep(self, outStar):
        self.subtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
        self.subtomoSet.copyInfo(self.inputSubtomos.get())
        # Read generated star file and create the output objects
        warningMsg, _ = readPysegSubtomograms(outStar,
                                              self.inputSubtomos.get(),
                                              self.subtomoSet)
        if warningMsg:
            self.warningMsg = String(warningMsg)
            self._store()

        self._defineOutputs(**{outputObjects.subtomograms.name: self.subtomoSet})
        self._defineSourceRelation(self.inputSubtomos.get(), self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('*Generated files location*:\n'
                           '\t- Subtomograms files directory: %s\n'
                           '\t- Star file: %s\n'
                           '%s' %
                           (self._getExtraPath(POST_REC_OUT),
                            self._getExtraPath(POST_REC_OUT + '.star'),
                            SEE_METHODS_TAB))
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            mbMask = self.mbMask.get()
            if mbMask:
                methods.append('*Membrane suppression applied with*\n'
                               '\t- MembraneMask = %s\n'
                               '\t- SuppressionFactor = %1.2f\n' %
                               (mbMask.getFileName(), self.mbSupFactor.get()))
            if self.doGaussLowPassFilter.get():
                methods.append('*Gaussian low pass filter applied with*\n'
                               '\t- CutOffResolution[nm] = %2.1f\n'
                               '\t- AmplitudeCutOff = %1.2f\n'
                               '\t- AppliedToCTF = %s\n' %
                               (self.cutOffRes.get(), self.ampCutOff.get(), self.filterCTF.get()))

        return methods

    def _validate(self):
        validationMsg = []
        tol = 0.01
        subTomosRes = self.inputSubtomos.get().getSamplingRate()
        inMask = self.inMask.get()
        maskRes = inMask.getSamplingRate()
        mbMask = self.mbMask.get()
        if abs(maskRes - subTomosRes) > tol:
            validationMsg.append('Sampling rate of the input subtomograms and the input mask should be the same\n'
                                 '%2.3f != %2.3f' % (subTomosRes, maskRes))
        if mbMask:
            mbMaskRes = mbMask.getSamplingRate()
            if abs(mbMaskRes - subTomosRes) > tol:
                validationMsg.append('Sampling rate of the input subtomograms and the input membrane suppression mask '
                                     'should be the same\n'
                                     '%2.3f != %2.3f' % (subTomosRes, mbMaskRes))
        return validationMsg

    # --------------------------- UTIL functions -----------------------------------

    def _getCommand(self, outDir, outStar):
        posRecCmd = ' '
        posRecCmd += '%s ' % Plugin.getHome(POST_REC_SCRIPT_MEMB_ATT)
        posRecCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        posRecCmd += '--inMask %s ' % getFinalMaskFileName(self.inMask.get())
        posRecCmd += '--inMaskMbSup %s ' % (getFinalMaskFileName(self.mbMask.get()) if self.mbMask.get() else 'None')
        posRecCmd += '--mbSupFactor %s ' % (self.mbSupFactor.get() if self.mbSupFactor.get() else '0')
        posRecCmd += '--doGaussLowPass %s ' % False
        posRecCmd += '--outDir %s ' % outDir
        posRecCmd += '--outStar %s ' % outStar
        posRecCmd += '-j %s ' % self.numberOfThreads.get()
        return posRecCmd

