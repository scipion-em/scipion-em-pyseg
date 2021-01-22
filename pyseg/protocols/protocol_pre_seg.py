from pwem.protocols import EMProtocol
from pyworkflow.protocol import FileParam, NumericListParam, IntParam, FloatParam, GT
from pyworkflow.utils import Message
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import PRESEG_SCRIPT


class ProtPySegPostRecParticles(EMProtocol, ProtTomoBase):
    """"""

    _label = 'Pre-process segmented CIRCULAR membranes'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inStar', FileParam,
                      label='Seg particles star file',
                      important=True,
                      allowsNull=False,
                      help='Star file obtained in PySeg seg step.')
        group = form.addGroup('Sub-volume splitting')
        group.addParam('spSplit', NumericListParam,
                       label='Number of splits (X, Y, Z)',
                       help='')  # TODO
        group.addParam('spOffVoxels', IntParam,
                       label='Offset voxels',
                       allowsNull=False,
                       validators=GT(0),
                       help='')  # TODO
        group = form.addGroup('Membrane segmentation')
        group.addParam('sgVoxelSize', FloatParam,
                       label='Voxel size (Å/voxel)',
                       validators=GT(0),
                       important=True,
                       allowsNull=False)
        group.addParam('sgThreshold', IntParam,
                       label='Density threshold',
                       validators=GT(0),
                       help='')  # TODO
        group.addParam('sgSizeThreshold', IntParam,
                       label='Size threshold (voxels)',
                       condition='sgThreshold > 0',
                       validators=GT(0),
                       help='')  # TODO
        group.addParam('sgMembThk', FloatParam,
                       label='Segmented membrane thickness (Å)',
                       allowsNull=False,
                       validators=GT(0))
        group.addParam('sgMembNeigh', FloatParam,
                       label='Segmented mebmrane neighbours (Å)',
                       allowsNull=False,
                       validators=GT(0),
                       help='')  # TODO

    def _insertAllSteps(self):
        self._insertFunctionStep('pysegPreSegStep')
        self._insertFunctionStep('getMembraneCenterStep')
        starWithCenters = self._insertFunctionStep('pysegPreSegCenteredStep')
        self._insertFunctionStep('createOutputStep', starWithCenters)

    def pysegPreSegStep(self):
        inStar = self.inStar.get()
        outDir = self._getTmpPath()

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPreSegCmd(inStar, outDir))

    def getMembraneCenterStep(self):
        pass

    def pysegPreSegCenteredStep(self, starWithCenters):
        outDir = self._getExtraPath()

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPreSegCmd(starWithCenters, outDir))

    def createOutputStep(self):
        pass
        # # Read generated star file and create the output objects
        # self.subtomoSet = self._createSetOfSubTomograms()
        # self.subtomoSet.setSamplingRate(self.inMask.get().getSamplingRate())
        # warningMsg = readStarFile(self, self.subtomoSet, RELION_SUBTOMO_STAR, starFile=outStar)
        # if warningMsg:
        #     self.warningMsg = String(warningMsg)
        #     self._store()
        #
        # self._defineOutputs(outputSubTomograms=self.subtomoSet)

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
        preSegCmd += '-spOffVoxels %s ' % self.spOffVoxels.get()
        preSegCmd += '-sgVoxelSize %s ' % self.sgVoxelSize.get()/10  # required in nm
        preSegCmd += '-sgThreshold %s ' % self.sgThreshold.get()
        preSegCmd += '-sgSizeThreshold %s ' % self.sgSizeThreshold.get()
        preSegCmd += '-sgMembThk %s ' % self.sgMembThk.get()/10  # required in nm
        preSegCmd += '-sgMembNeigh %s ' % self.sgMembThk.get()/10  # required in nm

        return preSegCmd
