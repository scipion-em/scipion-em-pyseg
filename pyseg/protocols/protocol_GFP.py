from os import environ
from os.path import join

from pwem.protocols import EMProtocol, FileParam
from pyworkflow.protocol import FloatParam, NumericListParam, EnumParam
from pyworkflow.utils import Message, makePath, removeBaseExt
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import GRAPHS_OUT, GRAPHS_SCRIPT, FILS_OUT, FILS_SCRIPT, FILS_SOURCES, FILS_TARGETS, \
    PICKING_OUT, PICKING_SCRIPT, PICKING_SLICES
from pyseg.convert import readStarFile


class ProtPySegGFP(EMProtocol, ProtTomoBase):
    """"""

    _label = 'Graphs-Fils-Picking'
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

        form.addSection(label='Graphs')
        group = form.addGroup('GraphMCF')
        group.addParam('sSig', FloatParam,
                       label='Sigma for gaussian filtering',
                       default=1,
                       allowsNull=False,
                       help='Sigma for Gaussian fltering input tomograms. It allows to smooth '
                            'small and irrelevant features and increases SNR.')
        group.addParam('vDen', FloatParam,
                       label='Vertex density within membranes (nm³)',
                       default=0.0035,
                       allowsNull=False,
                       help='Vertex density within membranes. It allows to adjust simplifcation '
                            'adaptively for every tomogram.')
        group.addParam('vRatio', FloatParam,
                       label='Avg ratio vertex/edge of graph within membrane',
                       default=4,
                       allowsNull=False,
                       help='Averaged ratio vertex/edge in the graph within membrane.')
        group.addParam('maxLen', FloatParam,
                       label='Shortest distance to membrane (nm)',
                       default=10,
                       allowsNull=False,
                       help='Maximum euclidean shortest distance to membrane in nm.')

        form.addSection(label='Fils')
        group = form.addGroup('Graph thresholding')
        group.addParam('thMode', EnumParam,
                       default=0,
                       choices=['in', 'out'],
                       label='Orientation with respect to the membrane/filament',
                       display=EnumParam.DISPLAY_HLIST)
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

        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        graphsStarFile = self._insertFunctionStep('pysegGraphs')
        filsStarFile = self._insertFunctionStep('pysegFils', graphsStarFile)
        pickingStarFile = self._insertFunctionStep('pysegPicking', filsStarFile)
        self._insertFunctionStep('createOutputStep')

    def pysegGraphs(self):
        # Generate output dir
        outDir = self._getExtraPath(GRAPHS_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, '%s %s %s %s %s %s %s %s %s' % (
            self._getPysegScript(GRAPHS_SCRIPT),
            self.inStar.get(),
            outDir,
            self.pixelSize.get(),
            self.numberOfThreads.get(),
            self.sSig.get(),  # Sigma for gaussian filtering
            self.vDen.get(),  # Vertex density within membranes (nm³)
            self.vRatio.get(),  # Avg ratio vertex/edge of graph within membrane
            self.maxLen.get()))  # Shortest distance to membrane (nm)

        # Generated star file
        return join(outDir, removeBaseExt(self.inStar.get()) + '_mb_graph.star')

    def pysegFils(self, graphsStarFile):
        # Generate output dir
        outDir = self._getExtraPath(FILS_OUT)
        makePath(outDir)

        # Manage the 2-element list inputs
        gRgEud = self._ListAsStr2ListOfNum(self.gRgEud.get())
        gRgLen = self._ListAsStr2ListOfNum(self.gRgLen.get())
        gRgSin = self._ListAsStr2ListOfNum(self.gRgSin.get())

        # Script called
        Plugin.runPySeg(self, PYTHON, '%s %s %s %s %s %s %s %s %s  %s %s %s' % (
            self._getPysegScript(FILS_SCRIPT),
            self._getExtraPath(GRAPHS_OUT, removeBaseExt(self.inStar.get()) + '_mb_graph.star'),
            outDir,
            Plugin.getHome(FILS_SOURCES),
            Plugin.getHome(FILS_TARGETS),
            self.thMode.get(),
            gRgEud[0],
            gRgEud[1],
            gRgLen[0],
            gRgLen[1],
            gRgSin[0],
            gRgSin[1]))

        # Generated star file
        return join(outDir, 'fil_' + removeBaseExt(FILS_SOURCES) + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')

    def pysegPicking(self, filsStarFile):
        # Generate output dir
        outDir = self._getExtraPath(PICKING_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, '%s %s %s %s %s %s' % (
            self._getPysegScript(PICKING_SCRIPT),
            self._getExtraPath(FILS_OUT, 'fil_' + removeBaseExt(FILS_SOURCES)
                               + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star'),
            outDir,
            Plugin.getHome(PICKING_SLICES),
            self.peakTh.get(),
            self.peakNs.get()))

        # return join(outDir, removeBaseExt(self.inStar.get()) + '_parts.star')

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
    @staticmethod
    def _getPysegScript(scriptName):
        return join(environ.get("SCIPION_HOME", None), Plugin.getHome(scriptName))

    @staticmethod
    def _ListAsStr2ListOfNum(inList):
        """Convert string of integers separated by spaces to a list of integers"""
        return [int(i) for i in inList.split()]