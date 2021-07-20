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

from pwem.protocols import EMProtocol, PointerParam
from pyworkflow import BETA
from pyworkflow.protocol import String, FloatParam, LEVEL_ADVANCED, BooleanParam, GT, LE
from pyworkflow.utils import Message, makePath
from reliontomo.convert import writeSetOfSubtomograms
from scipion.constants import PYTHON
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_OUT, POST_REC_SCRIPT_MEMB_ATT
from pyseg.convert import readStarFile, RELION_SUBTOMO_STAR


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
                       label='Mask for membrane supression')
        group.addParam('mbSupFactor', FloatParam,
                       label='Membrane suppression factor (0, 1]',
                       validators=[GT(0), LE(1)],
                       condition='mbMask')
        group = form.addGroup('Gaussian low pass filter', expertLevel=LEVEL_ADVANCED)
        group.addParam('doGaussLowPassFilter', BooleanParam,
                       label='Apply filter to the particles?',
                       default=False)
        group.addParam('cutOffRes', FloatParam,
                       label='Cut-off resolution (nm)',
                       condition='doGaussLowPassFilter')
        group.addParam('ampCutOff', FloatParam,
                       label='Amplitude at cut-off (0, 1]',
                       validators=[GT(0), LE(1)],
                       condition='doGaussLowPassFilter')
        group.addParam('filterCTF', BooleanParam,
                       label='Apply this filter to the CTF?',
                       default=False,
                       condition='doGaussLowPassFilter')
        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        outStar = self._getExtraPath(POST_REC_OUT + '.star')
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.pysegPostRec, outStar)
        self._insertFunctionStep(self.createOutputStep, outStar)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        imgSet = self.inputSubtomos.get()
        imgStar = self._getExtraPath(self.inStarName)
        writeSetOfSubtomograms(imgSet, imgStar, isPyseg=True)

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
        warningMsg = readStarFile(self, self.subtomoSet, RELION_SUBTOMO_STAR, starFile=outStar)
        if warningMsg:
            self.warningMsg = String(warningMsg)
            self._store()

        self._defineOutputs(outputSetOfSubtomogram=self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('Generated files location:\n'
                           'Subtomograms files directory: %s\n'
                           'Star file: %s' %
                           (self._getExtraPath(POST_REC_OUT),
                            self._getExtraPath(POST_REC_OUT + '.star')))
        return summary

    def _validate(self):
        validationMsg = []
        tol = 0.001
        subTomosRes = self.inputSubtomos.get().getSamplingRate()
        maskRes = self.inMask.get().GetSamplingRate()
        mbMaskRes = self.mbMask.get()
        if abs(maskRes - subTomosRes) > tol:
            validationMsg.append('Sampling rate of the input subtomograms and the input mask should be the same\n'
                                 '%2.3f != %2.3f' % (subTomosRes, maskRes))
        if mbMaskRes:
            if abs(mbMaskRes - subTomosRes) > tol:
                validationMsg.append('Sampling rate of the input subtomograms and the input membrane suppression mask '
                                     'should be the same\n'
                                     '%2.3f != %2.3f' % (subTomosRes, mbMaskRes))
        if self.doGaussLowPassFilter.get():
            failingInputs = []
            if not self.cutOffRes.get():
                failingInputs.append('Cut-off resolution')
            if not self.ampCutOff.get():
                failingInputs.append('Amplitude at cut-off')
            if failingInputs:
                validationMsg.append('Inputs *%s* should not be empty to apply the gaussian low pas filter' %
                                     ', '.join(failingInputs))
        return validationMsg

    # --------------------------- UTIL functions -----------------------------------

    def _getCommand(self, outDir, outStar):
        posRecCmd = ' '
        posRecCmd += '%s ' % Plugin.getHome(POST_REC_SCRIPT_MEMB_ATT)
        posRecCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        posRecCmd += '--inMask %s ' % self.inMask.get().getFileName()
        posRecCmd += '--inMaskMbSup %s ' % self.mbMask.get().getFileName()
        posRecCmd += '--mbSupFactor %s ' % self.mbSupFactor.get() if self.mbSupFactor.get() else None
        posRecCmd += '--doGaussLowPass %s ' % self.doGaussLowPassFilter.get()
        posRecCmd += '--resolution %s ' % self.inputSubtomos.get().getSamplingRate()
        posRecCmd += '--cutOffRes %s ' % self.cutOffRes.get()
        posRecCmd += '--ampCutOff %s ' % self.ampCutOff.get()
        posRecCmd += '--filterCTF %s ' % self.filterCTF.get()
        posRecCmd += '--outDir %s ' % outDir
        posRecCmd += '--outStar %s ' % outStar
        posRecCmd += '-j %s ' % self.numberOfThreads.get()
        return posRecCmd
