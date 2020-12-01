import json
from os import environ
from os.path import join

from pwem.protocols import EMProtocol, FileParam, PointerParam, StringParam
from pyworkflow.protocol import IntParam, GT, String, FloatParam, NumericListParam, PathParam
from pyworkflow.utils import Message, makePath
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin, BRANCH
from pyseg.constants import POST_REC_OUT
from pyseg.convert import readStarFile


class ProtPySegGFPRParticles(EMProtocol, ProtTomoBase):
    """"""

    _label = 'Graphs-Fils-Picking-Rec particles'
    # warningMsg = None
    # subtomoSet = None

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
        form.addParam('pixelSize', FloatParam,
                      label='Pixel size (nm/voxel)',
                      default=1,
                      important=True,
                      allowsNull=False,
                      help='Input tomograms voxel size (nm/voxel)')
        form.addParam('nMPI', IntParam,
                      pointerClass='VolumeMask',
                      label='Number of processors',
                      default=6,
                      validators=[GT(0)],
                      help='Multiprocessing settings, number of processors dedicated '
                           'to this protocol execution.')

        form.addSection(label='Graphs')
        form.addParam('sSig', FloatParam,
                      label='Sigma for gaussian filtering',
                      default=1,
                      allowsNull=False,
                      help='Sigma for Gaussian fltering input tomograms. It allows to smooth '
                           'small and irrelevant features and increases SNR.')
        form.addParam('vDen', FloatParam,
                      label='Vertex density within membranes (nm³)',
                      default=0.0035,
                      allowsNull=False,
                      help='Vertex density within membranes. It allows to adjust simplifcation '
                           'adaptively for every tomogram.')
        form.addParam('vRatio', FloatParam,
                      label='Avg ratio vertex/edge of graph within membrane',
                      default=4,
                      allowsNull=False,
                      help='Averaged ratio vertex/edge in the graph within membrane.')
        form.addParam('maxLen', FloatParam,
                      label='Shortest distance to membrane (nm)',
                      default=10,
                      allowsNull=False,
                      help='Maximum euclidean shortest distance to membrane in nm.')

        form.addSection(label='Fils')
        group = form.addGroup('Filament geometry')
        group.addParam('gRgEud', NumericListParam,
                       label='Euclidean distance of vertices source-target (nm)',
                       default='1 60',
                       allowsNull=False,
                       help='Euclidean distance between source and target vertices in nm.')
        group.addParam('gRgLen', NumericListParam,
                       label='Geodesic distance of vertices source-target (nm)',
                       default='1 25',
                       allowsNull=False,
                       help='Geodesic distance trough the graph between source and target vertices in nm.')
        group.addParam('gRgSin', NumericListParam,
                       label='Filament sinuosity',
                       default='0 3',
                       allowsNull=False,
                       help='Filament sinuosity, geodesic/euclidean distances ratio.')

        form.addSection(label='Picking')
        group = form.addGroup('Peaks cleaning')
        group.addParam('peakTh', FloatParam,
                       label='Percentile of points to discard by their density level.',
                       default=0,
                       allowsNull=False)
        group.addParam('peakNs', FloatParam,
                       label='Min distance beyween selected points (nm).',
                       default=0.5,
                       allowsNull=False,
                       help='Scale suppression in nm, two selected points cannot be closer than this distance.')

        form.addSection(label='Rec')
        form.addParam('missingWedgeCtf', PathParam,
                      label='Missing wedge CTF',
                      important=True,
                      allowsNull=False,
                      help='Missing wedge file path or directory (if more than 1 file is desired to be loaded).')
        form.addParam('filesPattern', StringParam,
                      label='Pattern',
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addParam('inMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask',
                      help='Mask used for the post processing. '
                           'Is is Required if field _psSegImage is not contained in'
                           'the star file generated in the picking step.')

    def _insertAllSteps(self):
        self._insertFunctionStep('pysegGraphs')
        self._insertFunctionStep('pysegFils')
        self._insertFunctionStep('pysegPicking')
        self._insertFunctionStep('pysegRec')
        self._insertFunctionStep('createOutputStep')

    def pysegGraphs(self):
        pass

    def pysegFilss(self):
        pass

    def pysegPicking(self):
        pass

    def pysegRec(self):
        pass

    def createOutputStep(self):
        pass
        # self._defineOutputs(outputSubTomograms=self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        pass
        # summary = []
        # if self.isFinished():
        #     summary.append('Generated files location:\n'
        #                    'Subtomograms files: *%s*\n'
        #                    'Star file: *%s*' %
        #                    (self._getExtraPath(POST_REC_OUT),
        #                     self._getExtraPath(POST_REC_OUT + '.star')))
        # return summary

    # --------------------------- UTIL functions -----------------------------------

