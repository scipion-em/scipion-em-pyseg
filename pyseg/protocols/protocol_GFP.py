from collections import OrderedDict
from os.path import basename, join, abspath
import xml.etree.ElementTree as ET

from pwem.protocols import EMProtocol
from pyworkflow.protocol import FloatParam, NumericListParam, EnumParam, PointerParam, IntParam
from pyworkflow.utils import Message, makePath, removeBaseExt, copyFile
from scipion.constants import PYTHON
from tomo.objects import SetOfCoordinates3D
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition

from pyseg import Plugin
from pyseg.constants import GRAPHS_OUT, GRAPHS_SCRIPT, FILS_OUT, FILS_SCRIPT, FILS_SOURCES, FILS_TARGETS, \
    PICKING_OUT, PICKING_SCRIPT, PICKING_SLICES
from pyseg.convert import readStarFile, PYSEG_PICKING_STAR, getTomoSetFromStar

TH_MODE_IN = 0
TH_MODE_OUT = 1

# Fils sources xml fields
SEG_LABEL_S = 'segLabelS'
MIN_EUC_DIST_S = 'minEucDistS'
MAX_EUC_DIST_S = 'maxEucDistS'
EUC_RANGE_S = 'eucRangeS'
MIN_GEO_DIST_S = 'minGeoDistS'
MAX_GEO_DIST_S = 'maxGeoDistS'
GEO_RANGE_S = 'geoDistRangeS'
MIN_GEO_LEN_S = 'minGeoLenS'
MAX_GEO_LEN_S = 'maxGeoLenS'
GEO_LEN_RANGE_S = 'geoLenRangeS'
MIN_FIL_SINU_S = 'minSinuS'
MAX_FIL_SINU_S = 'maxSinuS'
SINU_RANGE_S = 'sinurangeS'

# Fils targets xml fields
SEG_LABEL_T = 'segLabelT'
MIN_EUC_DIST_T = 'minEucDistT'
MAX_EUC_DIST_T = 'maxEucDistT'
EUC_RANGE_T = 'eucRangeT'
MIN_GEO_DIST_T = 'minGeoDistT'
MAX_GEO_DIST_T = 'maxGeoDistT'
GEO_RANGE_T = 'geoDistRangeT'
MIN_GEO_LEN_T = 'minGeoLenT'
MAX_GEO_LEN_T = 'maxGeoLenT'
GEO_LEN_RANGE_T = 'geoLenRangeT'
MIN_FIL_SINU_T = 'minSinuT'
MAX_FIL_SINU_T = 'maxSinuT'
SINU_RANGE_T = 'sinurangeT'


class ProtPySegGFP(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """"""

    _label = 'Graphs-Fils-Picking'
    _xmlSources = None
    _xmlTargets = None
    tomoSet = None
    warningMsg = None
    acquisitionParams = {
            'angleMin': 90,
            'angleMax': -90,
            'step': None,
            'angleAxis1': None,
            'angleAxis2': None
        }

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSegProt', PointerParam,
                      pointerClass='ProtPySegPreSegParticles',
                      label='Pre segmentation',
                      important=True,
                      allowsNull=False,
                      help='Pointer to preseg protocol.')
        # form.addParam('inStar', FileParam,
        #               label='Seg particles star file',
        #               important=True,
        #               allowsNull=False,
        #               help='Star file obtained in PySeg seg step.')
        form.addParam('pixelSize', FloatParam,
                      label='Pixel size (Å/voxel)',
                      default=1,
                      important=True,
                      allowsNull=False,
                      help='Input tomograms voxel size (Å/voxel)')

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

        form.addSection(label='Fils Sources')
        self._defineFilsXMLParams(form, self._getXMLSourcesDefaultVals())

        form.addSection(label='Fils Targets')
        self._defineFilsXMLParams(form, self._getXMLTargetsDefaultVals(), isSources=False)

        form.addSection(label='Fils Refinement')
        group = form.addGroup('Graph thresholding')
        group.addParam('thMode', EnumParam,
                       default=0,
                       choices=['in', 'out'],
                       label='Orientation with respect to the membrane/filament',
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('Filament geometry refinement')
        group.addParam('gRgEud', NumericListParam,
                       label='Euclidean distance of vertices source-target (nm)',
                       default='1 1000',
                       allowsNull=False,
                       help='Euclidean distance between source and target vertices in nm.')
        group.addParam('gRgLen', NumericListParam,
                       label='Geodesic distance of vertices source-target (nm)',
                       default='1 1000',
                       allowsNull=False,
                       help='Geodesic distance trough the graph between source and target vertices in nm.')
        group.addParam('gRgSin', NumericListParam,
                       label='Filament sinuosity',
                       default='0 1000',
                       allowsNull=False,
                       help='Filament sinuosity, geodesic/euclidean distances ratio.')

        form.addSection(label='Picking')
        group = form.addGroup('Peaks cleaning')
        group.addParam('peakTh', FloatParam,
                       label='Percentile of points to discard by their density level.',
                       default=0,
                       allowsNull=False)
        group.addParam('peakNs', FloatParam,
                       label='Min distance between selected points (nm).',
                       default=0.5,
                       allowsNull=False,
                       help='Scale suppression in nm, two selected points cannot be closer than this distance.')

        form.addParallelSection(threads=4, mpi=0)

    @staticmethod
    def _defineFilsXMLParams(form, d, isSources=True):
        """d is a disctionary with the default values"""
        sectionName = 'Sources - ' if isSources else 'Targets - '
        paramList = list(d.keys())
        valList = list(d.values())
        form.addParam(paramList[0], IntParam,
                      label='Segmentation label',
                      default=valList[0],
                      allowsNull=False,
                      help='Value of vertex property membrane segmentation (mb_seg).')
        group = form.addGroup('%sEuclidean distance to membrane (nm)' % sectionName)
        group.addParam(paramList[1], FloatParam,
                       label='Min',
                       default=valList[1],
                       allowsNull=False)
        group.addParam(paramList[2], FloatParam,
                       label='Max',
                       default=valList[2],
                       allowsNull=False)
        group.addParam(paramList[3], EnumParam,
                       default=valList[3],
                       choices=['[min, max]', '[-inf, min] U [max, +inf]'],
                       label='range',
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('%sGeodesic distance to membrane (nm)' % sectionName)
        group.addParam(paramList[4], FloatParam,
                       label='Min',
                       default=valList[4],
                       allowsNull=False)
        group.addParam(paramList[5], FloatParam,
                       label='Max',
                       default=valList[5],
                       allowsNull=False)
        group.addParam(paramList[6], EnumParam,
                       default=valList[6],
                       choices=['[min, max]', '[-inf, min] U [max, +inf]'],
                       label='range',
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('%sGeodesic length to membrane (nm)' % sectionName)
        group.addParam(paramList[7], FloatParam,
                       label='Min',
                       default=valList[7],
                       allowsNull=False)
        group.addParam(paramList[8], FloatParam,
                       label='Max',
                       default=valList[8],
                       allowsNull=False)
        group.addParam(paramList[9], EnumParam,
                       default=valList[9],
                       choices=['[min, max]', '[-inf, min] U [max, +inf]'],
                       label='range',
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('%sFilament sinuosity' % sectionName)
        group.addParam(paramList[10], FloatParam,
                       label='Min',
                       default=valList[10],
                       allowsNull=False)
        group.addParam(paramList[11], FloatParam,
                       label='Max',
                       default=valList[11],
                       allowsNull=False)
        group.addParam(paramList[12], EnumParam,
                       default=valList[12],
                       choices=['[min, max]', '[-inf, min] U [max, +inf]'],
                       label='range',
                       display=EnumParam.DISPLAY_HLIST)

    def _insertAllSteps(self):
        self._insertFunctionStep('pysegGraphs')
        self._insertFunctionStep('pysegFils')
        self._insertFunctionStep('pysegPicking')
        self._insertFunctionStep('createOutputStep')

    def pysegGraphs(self):
        # Generate output dir
        outDir = self._getExtraPath(GRAPHS_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getGraphsCommand(outDir))

    def pysegFils(self):
        # Generate output dir
        outDir = self._getExtraPath(FILS_OUT)
        makePath(outDir)

        # Generate sources xml
        self._createFilsXmlFile(Plugin.getHome(FILS_SOURCES), outDir)
        # Generate targets xml
        self._createFilsXmlFile(Plugin.getHome(FILS_TARGETS), outDir, isSource=False)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getFilsCommand(outDir))

    def pysegPicking(self):
        # Generate output dir
        outDir = self._getExtraPath(PICKING_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPickingCommand(outDir))

    def createOutputStep(self):
        pickingStarFile = self.getPickingStarFileName()
        samplingRate = self.pixelSize.get()

        tomoSet = self._createSetOfTomograms()
        tomoSet.setSamplingRate(samplingRate)
        self.tomoSet = tomoSet
        getTomoSetFromStar(self, pickingStarFile)
        self._defineOutputs(outputTomograms=self.tomoSet)

        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(tomoSet, suffix)
        coordsSet.setSamplingRate(samplingRate)
        # coordsSet.setPrecedents(tomoSet)
        readStarFile(self, coordsSet, PYSEG_PICKING_STAR, starFile=pickingStarFile, invert=True)

        self._defineOutputs(outputCoordinates=coordsSet)

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
    # @staticmethod
    # def _getPysegScript(scriptName):
    #     return Plugin.getHome(scriptName)
    #
    # @staticmethod
    # def _ListAsStr2ListOfNum(inList):
    #     """Convert string of integers separated by spaces to a list of integers"""
    #     return [int(i) for i in inList.split()]
    #
    def getPickingStarFileName(self):
        filsStar = self._getExtraPath(FILS_OUT, 'fil_' + removeBaseExt(FILS_SOURCES)
                                      + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')
        return self._getExtraPath(PICKING_OUT, removeBaseExt(filsStar) + '_parts.star')

    def _getGraphsCommand(self, outDir):
        pixSize = self.pixelSize.get()/10
        graphsCmd = ' '
        graphsCmd += '%s ' % Plugin.getHome(GRAPHS_SCRIPT)
        graphsCmd += '--inStar %s ' % self._getPreSegStarFile()  # self.inStar.get()
        graphsCmd += '--outDir %s ' % outDir
        graphsCmd += '--pixelSize %s ' % pixSize  # PySeg requires it in nm
        graphsCmd += '--sSig %s ' % self.sSig.get()
        graphsCmd += '--vDen %s ' % self.vDen.get()
        graphsCmd += '--veRatio %s ' % self.vRatio.get()
        graphsCmd += '--maxLen %s ' % self.maxLen.get()
        graphsCmd += '-j %s ' % self.numberOfThreads.get()
        return graphsCmd

    def _getFilsCommand(self, outDir):
        filsCmd = ' '
        filsCmd += '%s ' % Plugin.getHome(FILS_SCRIPT)
        filsCmd += '--inStar %s ' % self._getExtraPath(GRAPHS_OUT, removeBaseExt(self._getPreSegStarFile()) + '_mb_graph.star')
        filsCmd += '--outDir %s ' % outDir
        filsCmd += '--inSources %s ' % abspath(self._xmlSources)
        filsCmd += '--inTargets %s ' % abspath(self._xmlTargets)
        filsCmd += '--thMode %s ' % self._parseThModeSelection()
        filsCmd += '--gRgLen %s ' % self.gRgLen.get()
        filsCmd += '--gRgSin %s ' % self.gRgSin.get()
        filsCmd += '--gRgEud %s ' % self.gRgEud.get()
        return filsCmd

    def _getPickingCommand(self, outDir):
        pickingCmd = ' '
        pickingCmd += '%s ' % Plugin.getHome(PICKING_SCRIPT)
        pickingCmd += '--inStar %s ' % \
                      self._getExtraPath(FILS_OUT, 'fil_' + removeBaseExt(FILS_SOURCES) +
                                         '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')

        pickingCmd += '--outDir %s ' % outDir
        pickingCmd += '--slicesFile %s ' % Plugin.getHome(PICKING_SLICES)
        pickingCmd += '--peakTh %s ' % self.peakTh.get()
        pickingCmd += '--peakNs %s ' % self.peakNs.get()
        return pickingCmd

    def _parseThModeSelection(self):
        if self.thMode.get() == TH_MODE_IN:
            return 'in'
        else:
            return 'out'

    def _getPreSegStarFile(self):
        return self.inSegProt.get().getPresegOutputFile(self.inSegProt.get().getVesiclesCenteredStarFile())

    @staticmethod
    def _getXMLSourcesDefaultVals():
        d = OrderedDict()
        d[SEG_LABEL_S] = 1
        d[MIN_EUC_DIST_S] = 0
        d[MAX_EUC_DIST_S] = 15
        d[EUC_RANGE_S] = 0
        d[MIN_GEO_DIST_S] = 0
        d[MAX_GEO_DIST_S] = float('inf')
        d[GEO_RANGE_S] = 0
        d[MIN_GEO_LEN_S] = 0
        d[MAX_GEO_LEN_S] = 45
        d[GEO_LEN_RANGE_S] = 0
        d[MIN_FIL_SINU_S] = 0
        d[MAX_FIL_SINU_S] = 3
        d[SINU_RANGE_S] = 0
        return d

    @staticmethod
    def _getXMLTargetsDefaultVals():
        d = OrderedDict()
        d[SEG_LABEL_T] = 1
        d[MIN_EUC_DIST_T] = 0
        d[MAX_EUC_DIST_T] = 15
        d[EUC_RANGE_T] = 0
        d[MIN_GEO_DIST_T] = 0
        d[MAX_GEO_DIST_T] = float('inf')
        d[GEO_RANGE_T] = 0
        d[MIN_GEO_LEN_T] = 0
        d[MAX_GEO_LEN_T] = 45
        d[GEO_LEN_RANGE_T] = 0
        d[MIN_FIL_SINU_T] = 0
        d[MAX_FIL_SINU_T] = 3
        d[SINU_RANGE_T] = 0
        return d

    def _createFilsXmlFile(self, templateFile, outDir, isSource=True):
        EUCLIDEAN_DIST = 'eu_dst'
        GEODESIC_DIST = 'geo_dst'
        GEODESIC_LEN = 'geo_len'
        FIL_SINU = 'sin'

        # Copy the template xml to extra
        filename = join(outDir, basename(templateFile))
        copyFile(templateFile, filename)

        # Prepare the data to be read as expected by the xml file
        geomDict = OrderedDict()
        if isSource:
            self._xmlSources = filename
            segLabel = SEG_LABEL_S
            geomDict[EUCLIDEAN_DIST] = [MIN_EUC_DIST_S, MAX_EUC_DIST_S, EUC_RANGE_S]
            geomDict[GEODESIC_DIST] = [MIN_GEO_DIST_S, MAX_GEO_DIST_S, GEO_RANGE_S]
            geomDict[GEODESIC_LEN] = [MIN_GEO_LEN_S, MAX_GEO_LEN_S, GEO_LEN_RANGE_S]
            geomDict[FIL_SINU] = [MIN_FIL_SINU_S, MAX_FIL_SINU_S, SINU_RANGE_S]
        else:
            self._xmlTargets = filename
            segLabel = SEG_LABEL_T
            geomDict[EUCLIDEAN_DIST] = [MIN_EUC_DIST_T, MAX_EUC_DIST_T, EUC_RANGE_T]
            geomDict[GEODESIC_DIST] = [MIN_GEO_DIST_T, MAX_GEO_DIST_T, GEO_RANGE_T]
            geomDict[GEODESIC_LEN] = [MIN_GEO_LEN_T, MAX_GEO_LEN_T, GEO_LEN_RANGE_T]
            geomDict[FIL_SINU] = [MIN_FIL_SINU_T, MAX_FIL_SINU_T, SINU_RANGE_S]

        # Edit the corresponding fields
        xmlTree = ET.parse(filename)
        rootElement = xmlTree.getroot()
        mb_slice = rootElement.findall("mb_slice")
        mb_slice = mb_slice[0]
        mb_slice.find('side').text = str(getattr(self, segLabel).get())

        for key, valList in geomDict.items():
            for el in mb_slice.findall(key):
                if el.attrib['id'] == 'low':
                    el.text = str(getattr(self, valList[0]).get())
                elif el.attrib['id'] == 'high':
                    el.text = str(getattr(self, valList[1]).get())
                elif el.attrib['id'] == 'sign':
                    el.text = self._decodeRangeValue(getattr(self, valList[2]).get())

        # Write the modified xml file.
        xmlTree.write(filename, encoding='UTF-8', xml_declaration=True)

    @staticmethod
    def _decodeRangeValue(val):
        """Decode the range values and represent them as expected by pySeg"""
        # Choices are:
        #   0 --> [min, max], expected as '+'
        #   1 --> [-inf, min] U [max, +inf]'], expected as any other thing
        return '+' if val == 0 else '-'
