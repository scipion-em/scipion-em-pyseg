from pwem.protocols import EMProtocol, PointerParam
from pyworkflow import BETA
from pyworkflow.protocol import String
from pyworkflow.utils import Message, makePath
from reliontomo.convert import writeSetOfSubtomograms
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_SCRIPT, POST_REC_OUT
from pyseg.convert import readStarFile, RELION_SUBTOMO_STAR


class ProtPySegPostRecParticles(EMProtocol, ProtTomoBase):
    """"""

    _label = 'Post-process reconstructed particles'
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
        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        outStar = self._getExtraPath(POST_REC_OUT + '.star')
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('pysegPostRec', outStar)
        self._insertFunctionStep('createOutputStep', outStar)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        imgSet = self.inputSubtomos.get()
        imgStar = self._getExtraPath(self.inStarName)
        writeSetOfSubtomograms(imgSet, imgStar, isPysegPosRec=True)

    def pysegPostRec(self, outStar):
        # Generate output subtomo dir
        outDir = self._getExtraPath(POST_REC_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self. _getCommand(outDir, outStar))

    def createOutputStep(self, outStar):
        # Read generated star file and create the output objects
        self.subtomoSet = self._createSetOfSubTomograms()
        self.subtomoSet.setSamplingRate(self.inMask.get().getSamplingRate())
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

    # --------------------------- UTIL functions -----------------------------------

    def _getCommand(self, outDir, outStar):
        posRecCmd = ' '
        posRecCmd += '%s ' % Plugin.getHome(POST_REC_SCRIPT)
        posRecCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        posRecCmd += '--inMask %s ' % self.inMask.get().getFileName()
        posRecCmd += '--outDir %s ' % outDir
        posRecCmd += '--outStar %s ' % outStar
        posRecCmd += '-j %s ' % self.numberOfThreads.get()
        return posRecCmd
