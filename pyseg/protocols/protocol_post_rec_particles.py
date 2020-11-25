import json
from os import environ
from os.path import join

from pwem.protocols import EMProtocol, FileParam, PointerParam
from pyworkflow.protocol import IntParam, GT
from pyworkflow.utils import Message, makePath
from scipion.constants import PYTHON

from pyseg import Plugin, BRANCH


class ProtPySegPostRecParticles(EMProtocol):
    """"""

    _label = 'PySeg post-process reconstructed particles'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inStar', FileParam,
                      label='Particles star file',
                      important=True,
                      allowsNull=False,
                      help='Star file obtained in PySeg reconstruction.')
        form.addParam('inMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask',
                      important=True,
                      allowsNull=False,
                      help='Mask used for the post processing')
        form.addParam('nMPI', IntParam,
                      pointerClass='VolumeMask',
                      label='Number of processors',
                      default=6,
                      validators=[GT(0)],
                      help='Multiprocessing settings, number of processors dedicated '
                           'to this protocol execution.')

    def _insertAllSteps(self):
        self._insertFunctionStep('pysegPostRec')
        self._insertFunctionStep('createOutputStep')

    def pysegPostRec(self):
        outputName = 'subtomos_post_rec'
        makePath(outputName)
        # Inputs
        inStar = self.inStar.get()
        inMask = self.inMask.get().getFileName()
        outSubtomoDir = self._getExtraPath(outputName)
        outStar = self._getExtraPath(outputName + '.star')
        nMpi = self.nMPI.get()

        # Script called
        pyseg_post_rec = join(environ.get("SCIPION_HOME", None),
                              Plugin.getHome(join('pyseg_system-%s' % BRANCH, 'data', 'tutorials',
                                                  'synth_sumb', 'rln', 'post_rec_particles.py')))
        Plugin.runPySeg(self, PYTHON, '%s %s %s %s %s %s' % (pyseg_post_rec, inStar, inMask,
                                                             outSubtomoDir, outStar, nMpi))

    def createOutputStep(self):
        pass
        # train_data = CryocareTrainData(train_data=join(self.trainDataDir.get(), TRAIN_DATA_FN),
        #                                mean_std=join(self.trainDataDir.get(), MEAN_STD_FN),
        #                                patch_size=self._getPatchSize())
        # self._defineOutputs(train_data=train_data)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []

        # if self.isFinished():
        #     summary.append("Loaded training data info:\n"
        #                    "train_data_file = *{}*\n"
        #                    "normalization_file = *{}*\n"
        #                    "patch_size = *{}*".format(
        #                     join(self.trainDataDir.get(), TRAIN_DATA_FN),
        #                     join(self.trainDataDir.get(), MEAN_STD_FN),
        #                     self._getPatchSize()))
        return summary

    # --------------------------- UTIL functions -----------------------------------
    def _getPatchSize(self):
        with open(self.trainConfigFile.get()) as json_file:
            data = json.load(json_file)
            return data['patch_shape'][0]
