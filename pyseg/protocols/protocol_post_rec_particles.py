from pwem.protocols import EMProtocol, FileParam, PointerParam
from pyworkflow.protocol import String
from pyworkflow.utils import Message, makePath
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_SCRIPT, POST_REC_OUT
from pyseg.convert import readStarFile, RELION_SUBTOMO_STAR


class ProtPySegPostRecParticles(EMProtocol, ProtTomoBase):
    """"""

    _label = 'Post-process reconstructed particles'
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
        form.addParam('inputSubtomograms', PointerParam,
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
        self._insertFunctionStep('pysegPostRec', outStar)
        self._insertFunctionStep('createOutputStep', outStar)

    def pysegPostRec(self, outStar):
        # Generate output subtomo dir
        outDir = self._getExtraPath(POST_REC_OUT)
        makePath(outDir)
        # Script called
        pyseg_post_rec = Plugin.getHome(POST_REC_SCRIPT)
        Plugin.runPySeg(self, PYTHON, '%s %s %s %s %s %s' % (
            pyseg_post_rec,
            self.inStar.get(),  # In star file
            self.inMask.get().getFileName(),  # In mask
            outDir,  # Out subtomo dir
            outStar,  # Out star file
            self.numberOfThreads.get()))  # Number of MPI

    def createOutputStep(self, outStar):
        # Read generated star file and create the output objects
        self.subtomoSet = self._createSetOfSubTomograms()
        self.subtomoSet.setSamplingRate(self.inMask.get().getSamplingRate())
        warningMsg = readStarFile(self, self.subtomoSet, RELION_SUBTOMO_STAR, starFile=outStar)
        if warningMsg:
            self.warningMsg = String(warningMsg)
            self._store()

        self._defineOutputs(outputSubTomograms=self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('Generated files location:\n'
                           'Subtomograms files: *%s*\n'
                           'Star file: *%s*' %
                           (self._getExtraPath(POST_REC_OUT),
                            self._getExtraPath(POST_REC_OUT + '.star')))
        return summary

    # --------------------------- UTIL functions -----------------------------------

