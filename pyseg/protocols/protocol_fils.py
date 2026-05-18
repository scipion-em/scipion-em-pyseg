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
import glob
from collections import OrderedDict
from os import mkdir
from os.path import basename, join, abspath
import xml.etree.ElementTree as ET

from pwem.protocols import EMProtocol
from pyseg.convert.convert import splitPysegStarFile
from pyworkflow.protocol import FloatParam, NumericListParam, EnumParam, PointerParam, LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.utils import Message, removeBaseExt, copyFile, moveFile
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition

from pyseg import Plugin
from pyseg.constants import FILS_SCRIPT, FILS_SOURCES, FILS_TARGETS, MEMBRANE, \
    MEMBRANE_OUTER_SURROUNDINGS, PRESEG_AREAS_LIST, IN_STARS_DIR, OUT_STARS_DIR, FILS_OUT, GRAPHS_OUT, FILS_FILES
from pyseg.utils import encodePresegArea, createStarDirectories, genOutSplitStarFileName, getPrevPysegProtOutStarFiles

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


class ProtPySegFils(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """
    Extracts and refines filament networks from membrane-associated graph
    representations generated from segmented tomographic data. The protocol
    identifies biologically meaningful filament trajectories connecting
    selected membrane regions and filters them according to geometric and
    topological constraints relevant to cryo-electron tomography analysis.

    AI Generated:

    Filament Network Extraction (ProtPySegFils) - User Manual
        Overview

        The Filament Network Extraction protocol is designed to identify,
        refine, and characterize filamentous structures embedded within
        membrane-associated graph representations derived from segmented
        tomograms. Its primary objective is to reconstruct biologically
        meaningful connectivity patterns between membrane regions by tracing
        continuous filament paths through graph-based structural landscapes.

        In cryo-electron tomography workflows, this protocol is especially
        useful for studying cytoskeletal assemblies, membrane-associated
        scaffolds, extracellular filament systems, or protein networks that
        exhibit elongated and connected structural organization. Rather than
        focusing on isolated particles, the protocol emphasizes the recovery
        of continuous filament trajectories and their geometric properties.

        Biological Context and Applications

        Filamentous systems are central to many cellular processes,
        including intracellular transport, membrane shaping, structural
        reinforcement, signaling, and organelle organization. In many
        biological datasets, filaments appear as highly interconnected
        networks associated with membranes or membrane-adjacent regions.
        Understanding their organization often requires not only detecting
        individual densities but also reconstructing their continuity and
        spatial relationships.

        This protocol enables the user to define biologically relevant
        source and target regions from segmented membrane environments.
        Filaments are then evaluated according to spatial proximity,
        connectivity, curvature, and flexibility. This approach allows the
        protocol to isolate networks that are consistent with expected
        biological organization while suppressing noisy or structurally
        implausible paths.

        Inputs and Workflow

        The protocol operates on graph representations previously generated
        from segmented membrane data. These graphs contain spatial and
        topological information describing membrane-associated density
        landscapes within tomograms.

        The workflow begins by defining source and target membrane regions.
        These regions represent the biological compartments or membrane
        domains between which filament trajectories are expected to occur.
        Typical examples include membrane surfaces, lumenal regions,
        cytosolic regions, or extracellular surroundings.

        Once the source and target regions are selected, the protocol
        evaluates all possible graph trajectories connecting them and applies
        multiple geometric filters to retain only biologically meaningful
        filaments. The resulting filament network can then be used for
        particle picking, structural interpretation, connectivity analysis,
        or downstream quantitative studies.

        Source and Target Definition

        One of the most biologically important aspects of the protocol is
        the selection of source and target membrane regions. These choices
        define the biological context of the extracted network and strongly
        influence the resulting filament population.

        Source regions typically represent the membrane compartment from
        which filaments originate, while target regions define the
        destination environment. For example, a user may investigate
        trajectories extending from the membrane toward the cytosol or
        outward into extracellular regions.

        Proper biological interpretation depends heavily on choosing regions
        consistent with the expected organization of the specimen. Incorrect
        source-target combinations may generate filament populations lacking
        structural or functional relevance.

        Geometric Filtering and Refinement

        The protocol provides several geometric constraints that help refine
        filament selection and improve robustness. These constraints allow
        the user to focus on trajectories matching expected biological
        properties while excluding noisy graph connections.

        Euclidean distance constraints define the straight spatial proximity
        between filament vertices and membrane regions. These filters are
        useful for restricting filaments to biologically plausible distances
        from membranes or compartments.

        Geodesic distance and geodesic length constraints characterize the
        curved trajectories followed by filaments through the graph network.
        These parameters are particularly important when studying highly
        curved membrane systems or elongated filament assemblies.

        Filament sinuosity provides a measure of flexibility or tortuosity.
        Low sinuosity values correspond to relatively straight filaments,
        whereas higher values indicate more curved or flexible trajectories.
        From a biological perspective, this parameter can help distinguish
        rigid scaffold-like assemblies from flexible filament systems.

        Orientation and Membrane Relationship

        The protocol also allows the user to specify the orientation of
        filament analysis relative to the membrane environment. This becomes
        particularly important in systems where filament polarity or spatial
        directionality carries biological meaning.

        In membrane-associated systems, inward and outward orientations may
        correspond to distinct biological processes or structural
        compartments. Correct orientation selection helps ensure that the
        extracted network reflects the intended biological interpretation.

        Parallel Processing and Large Datasets

        Cryo-electron tomography datasets often contain large numbers of
        vesicles, membranes, or segmented regions. To improve scalability
        and computational efficiency, the protocol distributes processing
        across independent subsets of the data.

        This parallelized strategy is especially beneficial for large-scale
        studies involving dense filament systems or high-resolution tomograms.
        It allows complex datasets to be processed more efficiently while
        maintaining reproducibility and workflow continuity.

        Outputs and Interpretation

        The protocol produces refined filament network descriptions linking
        biologically selected membrane regions. These outputs can be used for
        downstream particle picking, structural interpretation, quantitative
        network analysis, or visualization within tomographic environments.

        Biologically, the extracted filaments should be interpreted as
        candidate connectivity pathways consistent with both the segmented
        membrane organization and the imposed geometric constraints. Their
        reliability depends strongly on segmentation quality, graph density,
        and the biological appropriateness of the selected filtering ranges.

        Practical Recommendations

        In most biological applications, it is advisable to begin with broad
        geometric ranges and progressively refine the parameters after visual
        inspection of the resulting filament network. Excessively restrictive
        thresholds may eliminate meaningful structures, whereas overly broad
        settings may retain noisy or biologically irrelevant trajectories.

        Euclidean and geodesic constraints are particularly important for
        balancing sensitivity and specificity. Flexible filament systems may
        require relaxed sinuosity constraints, while rigid structural
        assemblies often benefit from tighter geometric filtering.

        Careful validation against tomographic density maps is strongly
        recommended, especially when the extracted filaments will be used
        for quantitative biological interpretation or downstream structural
        analysis.

        Final Perspective

        Filament extraction from membrane-associated graph representations
        provides a powerful framework for studying structural connectivity
        within complex cellular environments. By combining membrane-aware
        region selection with geometric refinement criteria, the protocol
        enables biologically meaningful reconstruction of filament systems
        that would otherwise be difficult to isolate in noisy tomographic
        datasets.

        For most cryo-ET studies, successful filament analysis depends on
        thoughtful selection of membrane regions, biologically realistic
        geometric constraints, and careful interpretation of the resulting
        network organization within the structural context of the specimen.
    """

    _label = 'fils'
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        self._xmlSources = None
        self._xmlTargets = None
        self._inStarDir = None
        self._outStarDir = None

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inGraphsProt', PointerParam,
                      pointerClass='ProtPySegGraphs',
                      label='Graphs',
                      important=True,
                      allowsNull=False,
                      help='Pointer to graphs protocol.')

        form.addSection(label='Sources')
        self._defineFilsXMLParams(form, self._getXMLSourcesDefaultVals())

        form.addSection(label='Targets')
        self._defineFilsXMLParams(form, self._getXMLTargetsDefaultVals(), isSources=False)

        form.addSection(label='Refinement')
        group = form.addGroup('Graph thresholding')
        group.addParam('thMode', EnumParam,
                       default=TH_MODE_IN,
                       choices=['in', 'out'],
                       label='Orientation with respect to the membrane/filament',
                       display=EnumParam.DISPLAY_HLIST,
                       expertLevel=LEVEL_ADVANCED)
        group = form.addGroup('Filament geometry refinement ranges [min max]')
        group.addParam('gRgEud', NumericListParam,
                       label='Euclidean distance (STRAIGHT) range of vertices source-target (nm)',
                       default='1 1000',
                       allowsNull=False)
        group.addParam('gRgLen', NumericListParam,
                       label='Geodesic distance (CURVED) range of vertices source-target (nm)',
                       default='1 1000',
                       allowsNull=False)
        group.addParam('gRgSin', NumericListParam,
                       label='Filament sinuosity range (FLEXIBILITY, normally the ratio geoLen / eucLen)',
                       default='0 1000',
                       allowsNull=False)

        form.addParallelSection(threads=3, mpi=1)

    @staticmethod
    def _defineFilsXMLParams(form, d, isSources=True):
        """d is a disctionary with the default values"""
        sectionName = 'Sources - ' if isSources else 'Targets - '
        paramList = list(d.keys())
        valList = list(d.values())
        form.addParam(paramList[0], EnumParam,
                      choices=PRESEG_AREAS_LIST,
                      label='Filament area',
                      default=valList[0],
                      allowsNull=False,
                      help='Source or destination (depending if you are in the Sources or Targets tab) '
                           'area for the filament calculation.')
        group = form.addGroup('%sEuclidean (STRAIGHT) distance to membrane (nm)' % sectionName)
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
                       expertLevel=LEVEL_ADVANCED,
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('%sGeodesic distance to membrane (nm)' % sectionName, expertLevel=LEVEL_ADVANCED)
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
        group = form.addGroup('%sGeodesic (CURVED) length to membrane (nm)' % sectionName)
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
                       expertLevel=LEVEL_ADVANCED,
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('%sFilament sinuosity (FLEXIBILITY, normally the ratio '
                              'geodesicLen/euclideanLen)' % sectionName)
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
                       expertLevel=LEVEL_ADVANCED,
                       display=EnumParam.DISPLAY_HLIST)

    def _insertAllSteps(self):
        inStarDict = self._initialize()
        for starFile, outDir in inStarDict.items():
            self._insertFunctionStep(self.pysegFils, starFile, outDir, prerequisites=[],needsGPU=False)

    def _initialize(self):
        outDir = self._getExtraPath()
        # Split the input file into n (threads) files
        self._outStarDir, self._inStarDir = createStarDirectories(self._getExtraPath())
        # Generate sources xml
        self._createFilsXmlFile(Plugin.getHome(FILS_SOURCES), outDir)
        # Generate targets xml
        self._createFilsXmlFile(Plugin.getHome(FILS_TARGETS), outDir, isSource=False)
        # Generate 1 star file per vesicle to parallelize the calls to Fils and improve performance
        inStarFiles = []
        for inStar in sorted(glob.glob(self.inGraphsProt.get()._getExtraPath(OUT_STARS_DIR, '*.star'))):
            inStarFiles.extend(splitPysegStarFile(inStar, self._getExtraPath(IN_STARS_DIR),
                                                  j=1,
                                                  prefix=FILS_OUT + '_',
                                                  fileCounter=len(inStarFiles) + 1))
        # Associate a different output folder to each star file generated to store the fils resulting star file because
        # it is always generated with the same name, so there can be concurrency problems in parallelization
        inStarDict = {}
        filsResultsDir = self._getExtraPath(FILS_FILES)
        mkdir(filsResultsDir)
        for i, starFile in enumerate(inStarFiles):
            outDirName = join(filsResultsDir, 'outDir_%03d' % i)
            mkdir(outDirName)
            inStarDict[starFile] = outDirName

        return inStarDict

    def pysegFils(self, starFile, outDir):
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getFilsCommand(outDir, starFile))
        # Fils returns the same star file name, so it will be renamed to avoid overwriting
        moveFile(join(outDir, 'fil_mb_sources_to_no_mb_targets_net.star'),
                 genOutSplitStarFileName(self._outStarDir, starFile.replace(GRAPHS_OUT, FILS_OUT)))

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if self.isFinished():
            summary.append('*Filaments calculation*:\n\t- Source = %s\n\t- Target = %s\n' %
                           (PRESEG_AREAS_LIST[int(self.segLabelS.get())], PRESEG_AREAS_LIST[int(self.segLabelT.get())]))

        return summary

    # --------------------------- UTIL functions -----------------------------------
    def _getFilsCommand(self, outDir, starFile):
        filsCmd = ' '
        filsCmd += '%s ' % Plugin.getHome(FILS_SCRIPT)
        filsCmd += '--inStar %s ' % starFile
        filsCmd += '--outDir %s ' % outDir
        filsCmd += '--inSources %s ' % abspath(self._xmlSources)
        filsCmd += '--inTargets %s ' % abspath(self._xmlTargets)
        filsCmd += '--thMode %s ' % self._parseThModeSelection()
        filsCmd += '--gRgLen %s ' % self.gRgLen.get()
        filsCmd += '--gRgSin %s ' % self.gRgSin.get()
        filsCmd += '--gRgEud %s ' % self.gRgEud.get()
        return filsCmd

    def _getGraphsStarFile(self):
        prot = self.inGraphsProt.get()
        return prot._getExtraPath(removeBaseExt(prot._getPreSegStarFile()) + '_mb_graph.star')

    def _parseThModeSelection(self):
        if self.thMode.get() == TH_MODE_IN:
            return 'in'
        else:
            return 'out'

    @staticmethod
    def _getXMLSourcesDefaultVals():
        d = OrderedDict()
        d[SEG_LABEL_S] = MEMBRANE
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
        d[SEG_LABEL_T] = MEMBRANE_OUTER_SURROUNDINGS
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
            segLabel = encodePresegArea(self.segLabelS.get())
            geomDict[EUCLIDEAN_DIST] = [MIN_EUC_DIST_S, MAX_EUC_DIST_S, EUC_RANGE_S]
            geomDict[GEODESIC_DIST] = [MIN_GEO_DIST_S, MAX_GEO_DIST_S, GEO_RANGE_S]
            geomDict[GEODESIC_LEN] = [MIN_GEO_LEN_S, MAX_GEO_LEN_S, GEO_LEN_RANGE_S]
            geomDict[FIL_SINU] = [MIN_FIL_SINU_S, MAX_FIL_SINU_S, SINU_RANGE_S]
        else:
            self._xmlTargets = filename
            segLabel = encodePresegArea(self.segLabelT.get())
            geomDict[EUCLIDEAN_DIST] = [MIN_EUC_DIST_T, MAX_EUC_DIST_T, EUC_RANGE_T]
            geomDict[GEODESIC_DIST] = [MIN_GEO_DIST_T, MAX_GEO_DIST_T, GEO_RANGE_T]
            geomDict[GEODESIC_LEN] = [MIN_GEO_LEN_T, MAX_GEO_LEN_T, GEO_LEN_RANGE_T]
            geomDict[FIL_SINU] = [MIN_FIL_SINU_T, MAX_FIL_SINU_T, SINU_RANGE_S]

        # Edit the corresponding fields
        xmlTree = ET.parse(filename)
        rootElement = xmlTree.getroot()
        mb_slice = rootElement.findall("mb_slice")
        mb_slice = mb_slice[0]
        mb_slice.find('side').text = str(segLabel)

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

    # def _getGraphsOutStarFiles(prot):
    #     inStarList = glob.glob(prot.inGraphsProt.get()._getExtraPath(join(OUT_STARS_DIR, '*.star')))
    #     outStarFiles = []
    #     for inStarFile in inStarList:
    #         outStarFile = prot._getExtraPath(IN_STARS_DIR, basename(inStarFile))
    #         symlink(abspath(inStarFile), abspath(outStarFile))
    #         outStarFiles.append(outStarFile)
    #
    #     return outStarFiles
