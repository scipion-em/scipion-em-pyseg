from collections import OrderedDict
from os.path import basename, join, abspath
import xml.etree.ElementTree as ET

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FloatParam, NumericListParam, EnumParam, PointerParam, IntParam, FileParam, StringParam
from pyworkflow.utils import Message, makePath, removeBaseExt, copyFile
from scipion.constants import PYTHON
from tomo.objects import SetOfCoordinates3D
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition

from pyseg import Plugin
from pyseg.constants import GRAPHS_OUT, GRAPHS_SCRIPT, FILS_OUT, FILS_SCRIPT, FILS_SOURCES, FILS_TARGETS, \
    PICKING_OUT, PICKING_SCRIPT, PICKING_SLICES, FROM_SCIPION, FROM_STAR_FILE
from pyseg.convert import readStarFile, PYSEG_PICKING_STAR, getTomoSetFromStar

TH_MODE_IN = 0
TH_MODE_OUT = 1

# Fils slices xml fields
SIDE = 'side'
CONT = 'cont'

# CONT param codification
CUTTING_POINT = 0
PROJECTIONS = 1


class ProtPySegPicking(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """"""

    _label = 'Picking'
    _devStatus = BETA
    tomoSet = None
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
        form.addParam('filsFrom', EnumParam,
                      choices=['Scipion Protocol', 'Star file'],
                      default=0,
                      label='Choose graphs data source',
                      important=True,
                      display=EnumParam.DISPLAY_HLIST)
        form.addParam('inFilsProt', PointerParam,
                      pointerClass='ProtPySegFils',
                      label='Pre segmentation',
                      condition='filsFrom == 0',
                      important=True,
                      allowsNull=False,
                      help='Pointer to fils protocol.')
        form.addParam('inStar', FileParam,
                      label='Seg particles star file',
                      condition='filsFrom == 1',
                      important=True,
                      allowsNull=False,
                      help='Star file obtained in PySeg seg step.')
        form.addParam('pixelSize', FloatParam,
                      label='Pixel size (Å/voxel)',
                      default=1,
                      important=True,
                      allowsNull=False,
                      help='Input tomograms voxel size (Å/voxel)')

        form.addSection(label='Picking')
        self._defineFilsXMLParams(form, self._getSlicesXMLDefaultVals())

        form.addSection(label='Refinement')
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

    @staticmethod
    def _defineFilsXMLParams(form, d):
        """d is a disctionary with the default values"""
        paramList = list(d.keys())
        valList = list(d.values())
        form.addParam(paramList[0], IntParam,
                      label='Segmentation index area for picking',
                      default=valList[0],
                      allowsNull=False)
        form.addParam(paramList[1], EnumParam,
                      choices=['Cutting point', 'Cutting point + projections'],
                      default=0,
                      label='Find on two surfaces',
                      display=EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.")

    @staticmethod
    def _getSlicesXMLDefaultVals():
        d = OrderedDict()
        d[SIDE] = 1
        d[CONT] = CUTTING_POINT
        return d

    def _insertAllSteps(self):
        self._insertFunctionStep('pysegPicking')
        self._insertFunctionStep('createOutputStep')

    def pysegPicking(self):
        # Generate output dir
        outDir = self._getExtraPath()
        makePath(outDir)

        # Generate slices xml
        self._createPickingXmlFile(Plugin.getHome(PICKING_SLICES), outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPickingCommand(outDir))

    def createOutputStep(self):
        pickingStarFile = self._getPickingStarFileName()
        samplingRate = self.pixelSize.get()

        tomoSet = self._createSetOfTomograms()
        tomoSet.setSamplingRate(samplingRate)
        self.tomoSet = tomoSet
        getTomoSetFromStar(self, pickingStarFile)
        self._defineOutputs(outputTomograms=self.tomoSet)

        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(tomoSet, suffix)
        coordsSet.setSamplingRate(samplingRate)
        readStarFile(self, coordsSet, PYSEG_PICKING_STAR, starFile=pickingStarFile, invert=True)

        self._defineOutputs(outputCoordinates=coordsSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        pass

    # --------------------------- UTIL functions -----------------------------------
    def _getFilsStarFileName(self):
        source = self.filsFrom.get()
        if source == FROM_SCIPION:
            prot = self.inFilsProt.get()
            return prot._getExtraPath('fil_' + removeBaseExt(FILS_SOURCES)
                                      + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')

        elif source == FROM_STAR_FILE:
            return self.inStar.get()

    def _getPickingStarFileName(self):
        filsStar = self._getExtraPath('fil_' + removeBaseExt(FILS_SOURCES)
                                      + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')
        return self._getExtraPath(removeBaseExt(filsStar) + '_parts.star')

    def _getPickingCommand(self, outDir):
        pickingCmd = ' '
        pickingCmd += '%s ' % Plugin.getHome(PICKING_SCRIPT)
        pickingCmd += '--inStar %s ' % self._getFilsStarFileName()
        pickingCmd += '--outDir %s ' % outDir
        pickingCmd += '--slicesFile %s ' % self._xmlSlices
        pickingCmd += '--peakTh %s ' % self.peakTh.get()
        pickingCmd += '--peakNs %s ' % self.peakNs.get()
        return pickingCmd

    def _createPickingXmlFile(self, templateFile, outDir):
        # Copy the template xml to extra
        filename = join(outDir, basename(templateFile))
        self._xmlSlices = filename
        copyFile(templateFile, filename)

        # Edit the corresponding fields
        xmlTree = ET.parse(filename)
        rootElement = xmlTree.getroot()
        mb_slice = rootElement.findall("mb_slice")
        mb_slice = mb_slice[0]
        mb_slice.find(SIDE).text = str(getattr(self, SIDE).get())
        mb_slice.find(CONT).text = self._decodeContValue(getattr(self, CONT).get())

        # Write the modified xml file.
        xmlTree.write(filename, encoding='UTF-8', xml_declaration=True)

    @staticmethod
    def _decodeContValue(val):
        """Decode the cont values and represent them as expected by pySeg"""
        return '+' if val == CUTTING_POINT else '-'
