import json
from os import environ
from os.path import join

from pwem.protocols import EMProtocol, FileParam, PointerParam
from pyworkflow.protocol import IntParam, GT, String
from pyworkflow.utils import Message, makePath
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin, BRANCH
from pyseg.constants import POST_REC_OUT
from pyseg.convert import readStarFile


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
        # JORGE
        import os
        fname = "/home/jjimenez/Desktop/test_JJ.txt"
        if os.path.exists(fname):
            os.remove(fname)
        fjj = open(fname, "a+")
        fjj.write('JORGE--------->onDebugMode PID {}'.format(os.getpid()))
        fjj.close()
        print('JORGE--------->onDebugMode PID {}'.format(os.getpid()))
        import time
        time.sleep(10)
        # JORGE_END
        self._insertFunctionStep('pysegPostRec')
        self._insertFunctionStep('createOutputStep')

    def pysegPostRec(self):
        makePath(POST_REC_OUT)
        # Inputs
        inStar = self.inStar.get()
        inMask = self.inMask.get().getFileName()
        outSubtomoDir = self._getExtraPath(POST_REC_OUT)
        outStar = self._getExtraPath(POST_REC_OUT + '.star')
        nMpi = self.nMPI.get()

        # Script called
        pyseg_post_rec = join(environ.get("SCIPION_HOME", None),
                              Plugin.getHome(join('pyseg_system-%s' % BRANCH, 'data', 'tutorials',
                                                  'synth_sumb', 'rln', 'post_rec_particles.py')))
        Plugin.runPySeg(self, PYTHON, '%s %s %s %s %s %s' % (pyseg_post_rec, inStar, inMask,
                                                             outSubtomoDir, outStar, nMpi))

        # Read generated star file and create the output objects
        self.subtomoSet = self._createSetOfSubTomograms()
        self.subtomoSet.setSamplingRate(self.inMask.get().getSamplingRate())
        warningMsg = readStarFile(self, self.subtomoSet, starFile=outStar)
        if warningMsg:
            self.warningMsg = String(warningMsg)
            self._store()

    def createOutputStep(self):
        self._defineOutputs(outputSubTomograms=self.subtomoSet)

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

