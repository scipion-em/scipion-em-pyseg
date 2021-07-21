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
from os.path import exists, abspath

from pwem.protocols import EMProtocol, PointerParam
from pyworkflow import BETA
from pyworkflow.protocol import String, FloatParam, LEVEL_ADVANCED, BooleanParam, GT, LE, GE
from pyworkflow.utils import Message, makePath
from reliontomo.convert import writeSetOfSubtomograms
from scipion.constants import PYTHON
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_OUT, POST_REC_SCRIPT_MEMB_ATT, SEE_METHODS_TAB
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
        group = form.addGroup('Gaussian low pass filter', expertLevel=LEVEL_ADVANCED)
        group.addParam('doGaussLowPassFilter', BooleanParam,
                       label='Apply filter to the particles?',
                       default=False)
        group.addParam('cutOffRes', FloatParam,
                       label='Cut-off resolution (nm)',
                       condition='doGaussLowPassFilter')
        group.addParam('ampCutOff', FloatParam,
                       label='Amplitude at cut-off (0, 1]',
                       default=0.01,
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
        tol = 0.001
        subTomosRes = self.inputSubtomos.get().getSamplingRate()
        maskRes = self.inMask.get().getSamplingRate()
        mbMaskRes = self.mbMask.get()
        if abs(maskRes - subTomosRes) > tol:
            validationMsg.append('Sampling rate of the input subtomograms and the input mask should be the same\n'
                                 '%2.3f != %2.3f' % (subTomosRes, maskRes))
        if mbMaskRes:
            if abs(mbMaskRes.getSamplingRate() - subTomosRes) > tol:
                validationMsg.append('Sampling rate of the input subtomograms and the input membrane suppression mask '
                                     'should be the same\n'
                                     '%2.3f != %2.3f' % (subTomosRes, mbMaskRes))
        if self.doGaussLowPassFilter.get():
            # CTF3D required
            COORD_CTF_3D = '_3dcftMrcFile'
            coord3D = self.inputSubtomos.get().getCoordinates3D().get().getFirstItem()
            ctf3d = getattr(coord3D, COORD_CTF_3D, None)
            if ctf3d:
                ctf3d = abspath(ctf3d.get())
                if not exists(ctf3d):
                    validationMsg.append('CTF3D not found.\nValue read for the first coordinate is\n%s = %s'
                                         % (COORD_CTF_3D, ctf3d))
            else:
                validationMsg.append('CTF3D is required to apply the gaussian low pass filter. Please calculate it\n'
                                     'or do not apply the filter.')
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
        doGaussianLPFilter = self.doGaussLowPassFilter.get()
        cutOffRes = 0
        ampCutOff = 0
        if doGaussianLPFilter:
            cutOffRes = self.cutOffRes.get()
            ampCutOff = self.ampCutOff.get()
        posRecCmd = ' '
        posRecCmd += '%s ' % Plugin.getHome(POST_REC_SCRIPT_MEMB_ATT)
        posRecCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        posRecCmd += '--inMask %s ' % self.inMask.get().getFileName()
        posRecCmd += '--inMaskMbSup %s ' % (self.mbMask.get().getFileName() if self.mbMask.get() else 'None')
        posRecCmd += '--mbSupFactor %s ' % (self.mbSupFactor.get() if self.mbSupFactor.get() else '0')
        posRecCmd += '--doGaussLowPass %s ' % doGaussianLPFilter
        posRecCmd += '--resolution %s ' % (float(self.inputSubtomos.get().getSamplingRate()) / 10)  # in nm
        posRecCmd += '--cutOffRes %s ' % cutOffRes
        posRecCmd += '--ampCutOff %s ' % ampCutOff
        posRecCmd += '--filterCTF %s ' % self.filterCTF.get()
        posRecCmd += '--outDir %s ' % outDir
        posRecCmd += '--outStar %s ' % outStar
        posRecCmd += '-j %s ' % self.numberOfThreads.get()
        return posRecCmd
