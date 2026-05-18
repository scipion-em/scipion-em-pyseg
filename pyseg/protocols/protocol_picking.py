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
from collections import OrderedDict
from enum import Enum
from os.path import basename, join
import xml.etree.ElementTree as ET
from pwem.protocols import EMProtocol
from pyseg.convert import readPysegCoordinates
from pyworkflow.protocol import FloatParam, EnumParam, PointerParam, IntParam, LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.utils import Message, removeBaseExt, copyFile, moveFile
from scipion.constants import PYTHON
from tomo.objects import SetOfCoordinates3D, SetOfTomograms
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition
from pyseg import Plugin
from pyseg.constants import FILS_SOURCES, FILS_TARGETS, PICKING_SCRIPT, PICKING_SLICES, PRESEG_AREAS_LIST, MEMBRANE, \
    OUT_STARS_DIR, IN_STARS_DIR, FILS_OUT, PICKING_OUT

# Fils slices xml fields
from pyseg.utils import encodePresegArea, getPrevPysegProtOutStarFiles, createStarDirectories
from tomo.utils import getObjFromRelation

SIDE = 'side'
CONT = 'cont'

# CONT param codification
CUTTING_POINT = 0
PROJECTIONS = 1


class outputObjects(Enum):
    coordinates = SetOfCoordinates3D


class ProtPySegPicking(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """
    Extracts particle coordinates from filament networks associated with
    oriented membrane graphs in cryo-electron tomography datasets. The
    protocol identifies biologically relevant positions along filamentous
    structures and generates coordinate sets suitable for subtomogram
    extraction and downstream structural analysis.

    AI Generated:

    PySeg Picking (ProtPySegPicking) — User Manual
        Overview

        The PySeg Picking protocol is designed to identify and extract
        particle coordinates from filament networks embedded within or
        associated with segmented membrane environments. In cryo-electron
        tomography workflows, filamentous assemblies frequently organize
        along vesicles, organelles, or membrane surfaces, making targeted
        coordinate detection essential for subsequent subtomogram
        averaging and structural interpretation.

        This protocol is particularly useful for studies involving
        cytoskeletal filaments, membrane-associated polymers, tubular
        scaffolds, or molecular assemblies distributed along segmented
        membrane geometries. Rather than performing generic particle
        detection, the protocol focuses on biologically meaningful
        regions derived from previously analyzed filament networks and
        membrane graph representations.

        Inputs and Biological Context

        The protocol requires a previously generated filament network
        analysis as its main input. These filament networks define the
        structural organization used for coordinate detection and provide
        the geometric framework from which particles are selected. The
        workflow assumes that segmentation and graph construction have
        already been completed successfully.

        Optionally, users may specify a tomogram set to define the
        coordinate reference system. When no explicit tomogram set is
        provided, the protocol automatically associates coordinates with
        the tomograms linked to the original membrane segmentation
        workflow. This behavior simplifies integration within larger
        cryo-ET processing pipelines.

        The selected box size determines the extraction region that will
        later surround each picked coordinate. In biological practice,
        the box should be large enough to contain the complete molecular
        complex and nearby contextual information while avoiding
        excessive solvent or unrelated structures.

        Picking Strategy and Membrane Context

        A key feature of the protocol is the ability to define the
        membrane-related region where coordinate detection is performed.
        Depending on the biological question, users may focus on the
        membrane itself, nearby regions, or specific directional
        projections associated with filament geometry.

        This flexibility is important because many filament systems are
        asymmetrically distributed relative to membranes. For example,
        some molecular complexes bind directly to the membrane surface,
        whereas others extend outward into the cytosol or lumen. By
        controlling the picking region, the protocol enables biologically
        focused particle extraction adapted to the spatial organization
        of the sample.

        The protocol also supports differentiation between opposite sides
        of the membrane environment. This is particularly useful in
        systems where orientation matters biologically, such as polarized
        assemblies, membrane remodeling events, or asymmetric filament
        organization.

        Coordinate Refinement and Cleaning

        After candidate positions are identified, the protocol applies
        refinement procedures intended to improve the quality and
        interpretability of the resulting coordinate set. Density-based
        cleaning removes weak or unreliable detections, helping reduce
        false positives and improve downstream averaging consistency.

        A spatial suppression criterion is also available to avoid
        selecting coordinates that are unrealistically close to one
        another. In crowded filament networks, nearby detections may
        correspond to redundant representations of the same biological
        feature. Enforcing a minimal separation distance improves the
        statistical independence of extracted particles and reduces
        oversampling artifacts.

        From a biological perspective, these refinement operations are
        especially valuable when studying dense molecular environments
        where structural repetition and membrane curvature can complicate
        automatic detection.

        Outputs and Interpretation

        The protocol produces a set of three-dimensional coordinates
        associated with the selected tomograms. These coordinates define
        the spatial locations of candidate particles suitable for
        subtomogram extraction, averaging, classification, or structural
        mapping workflows.

        The generated coordinate set preserves the tomographic reference
        frame and remains linked to the biological context of the
        original segmentation and filament analysis. This allows users
        to correlate picked particles with membrane geometry, filament
        orientation, and spatial organization within the tomogram.

        In practical cryo-ET studies, the resulting coordinates often
        serve as the starting point for iterative structural refinement
        workflows. Visual inspection of representative coordinates is
        strongly recommended to verify biological plausibility and ensure
        that particle selection matches the intended structural targets.

        Practical Recommendations

        For most biological applications, users should begin with
        conservative cleaning thresholds and visually inspect the
        resulting coordinates before increasing filtering stringency.
        Excessive suppression may remove biologically meaningful
        particles, particularly in heterogeneous or low-contrast
        datasets.

        The selected membrane region should reflect the known biological
        organization of the system under study. Membrane-associated
        complexes with strong polarity or directional organization often
        benefit from side-specific picking strategies.

        Box size optimization is also important. Small boxes may truncate
        extended complexes or flexible regions, while excessively large
        boxes increase background variability and computational cost in
        downstream processing.

        Final Perspective

        In cryo-electron tomography workflows, coordinate picking is not
        simply a geometric detection task but a biologically informed
        selection process that strongly influences downstream structural
        interpretation. Careful definition of membrane context, filament
        organization, and coordinate refinement criteria enables more
        reliable particle extraction and improves the quality of
        subtomogram averaging and classification analyses.
    """

    _label = 'picking'

    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        self._tomoSet = None
        self.stepsExecutionMode = STEPS_PARALLEL
        self.tomoSet = None
        self.acquisitionParams = {
                'angleMin': 90,
                'angleMax': -90,
                'step': None
        }
        self._inStarDir = None
        self._outStarDir = None
        self._xmlSlices = None
        self._outStarFilesList = []

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inFilsProt', PointerParam,
                      pointerClass='ProtPySegFils',
                      label='Fils protocol',
                      important=True,
                      allowsNull=False,
                      help='Pointer to fils protocol.')
        form.addParam('inTomoSet', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms to refer the coordinates',
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='Tomograms to which the coordinates should be referred to. If empty, it is assumed that '
                           'the tomograms are the same from were the vesicles were segmented (pre-seg).')
        form.addParam('boxSize', IntParam,
                      label='Box size (pixels)',
                      default=20,
                      important=True,
                      allowsNull=False,
                      expertLevel=LEVEL_ADVANCED)

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

        form.addParallelSection(threads=3, mpi=1)

    @staticmethod
    def _defineFilsXMLParams(form, d):
        """d is a disctionary with the default values"""
        paramList = list(d.keys())
        valList = list(d.values())
        form.addParam(paramList[0], EnumParam,
                      label='Segmentation area for picking',
                      choices=PRESEG_AREAS_LIST,
                      default=valList[0],
                      allowsNull=False,
                      help='Area in which the cutting point or cutting point + projections of the '
                           'filament will be considered for the picking coordinates.')
        form.addParam(paramList[1], EnumParam,
                      choices=['Cutting point', 'Projected local minima'],
                      default=0,
                      label='Find on two surfaces',
                      display=EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.")

    @staticmethod
    def _getSlicesXMLDefaultVals():
        d = OrderedDict()
        d[SIDE] = MEMBRANE
        d[CONT] = CUTTING_POINT
        return d

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        allOutputId = []
        starFileList = self._convertInputStep()
        for starFile in starFileList:
            pId = self._insertFunctionStep(self.pysegPicking, starFile, prerequisites=[])
            allOutputId.append(pId)
        self._insertFunctionStep(self.createOutputStep, prerequisites=allOutputId)

    def _convertInputStep(self):
        outDir = self._getExtraPath()
        # Generate directories for input and output star files
        self._outStarDir, self._inStarDir = createStarDirectories(outDir)
        # Generate slices xml
        self._createPickingXmlFile(Plugin.getHome(PICKING_SLICES), outDir)
        # Get the star files generated in the fils protocol
        return getPrevPysegProtOutStarFiles(self.inFilsProt.get()._getExtraPath(OUT_STARS_DIR),
                                            self._getExtraPath(IN_STARS_DIR))

    def pysegPicking(self, starFile):
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPickingCommand(starFile))
        # Move output files to the corresponding directory
        outFile = self._getExtraPath(removeBaseExt(starFile) + '_parts.star')
        newFileName = join(self._outStarDir, basename(outFile).replace(FILS_OUT, PICKING_OUT))
        self._outStarFilesList.append(newFileName)
        moveFile(outFile, newFileName)

    def createOutputStep(self):
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(self._getTomoSet(), suffix)
        coordsSet.setSamplingRate(self._getSamplingRate())
        coordsSet.setBoxSize(self.boxSize.get())

        # Read the data from all the out star files
        tomoSet = self._getTomoSet()
        for outStar in self._outStarFilesList:
            readPysegCoordinates(outStar, coordsSet, tomoSet)

        if not coordsSet:
            raise Exception('ERROR! No coordinates were picked.')
        self._defineOutputs(**{outputObjects.coordinates.name: coordsSet})
        self._defineSourceRelation(self._getTomoSet(), coordsSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        valMsg = []
        if not self._getTomoFromRelations():
            valMsg.append("Unable to find the corresponding tomograms using the relations.")
        return valMsg

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if self.isFinished():
            outputCoords = getattr(self, outputObjects.coordinates.name, None)
            summary.append('*Picking*:\n\t- Picking area = %s\n\t- Particles picked = %i\n' %
                           (PRESEG_AREAS_LIST[int(self.side.get())], outputCoords.getSize()))

        return summary

    # --------------------------- UTIL functions -----------------------------------
    def _getFilsStarFileName(self):
        prot = self.inFilsProt.get()
        return prot._getExtraPath('fil_' + removeBaseExt(FILS_SOURCES)
                                  + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')

    def _getPickingStarFileName(self):
        filsStar = self._getExtraPath('fil_' + removeBaseExt(FILS_SOURCES)
                                      + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')
        return self._getExtraPath(removeBaseExt(filsStar) + '_parts.star')

    def _getPickingCommand(self, starFile):
        pickingCmd = ' '
        pickingCmd += '%s ' % Plugin.getHome(PICKING_SCRIPT)
        pickingCmd += '--inStar %s ' % starFile
        pickingCmd += '--outDir %s ' % self._getExtraPath()
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
        mb_slice.find(SIDE).text = str(encodePresegArea(getattr(self, SIDE).get()))
        mb_slice.find(CONT).text = self._decodeContValue(getattr(self, CONT).get())

        # Write the modified xml file.
        xmlTree.write(filename, encoding='UTF-8', xml_declaration=True)

    @staticmethod
    def _decodeContValue(val):
        """Decode the cont values and represent them as expected by pySeg"""
        return '+' if val == CUTTING_POINT else '-'

    def _getSamplingRate(self):
        return self._getTomoSet().getFirstItem().getSamplingRate()

    def _getTomoSet(self):
        if self._tomoSet is None:
            if self.inTomoSet.get():
                self._tomoSet = self.inTomoSet.get()
            else:
                self._tomoSet = self._getTomoFromRelations()
        return self._tomoSet

    def _getTomoFromRelations(self):
        # Get the tomograms climbing from this point of the workflow until the pre-seg and if there aren't tomograms
        # at that point, use the relations to go to the corresponding tomograms
        presegProt = self.inFilsProt.get().inGraphsProt.get().inSegProt.get()
        tomoSet = getObjFromRelation(presegProt.inTomoMasks.get(), self, SetOfTomograms)
        return tomoSet
