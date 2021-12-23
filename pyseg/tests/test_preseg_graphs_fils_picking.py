from glob import glob
from os import remove
from os.path import join, abspath, exists

from imod.protocols import ProtImodTomoNormalization
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from relion.convert import Table
from pyseg.constants import TOMOGRAM, VESICLE, PYSEG_LABEL, SEGMENTATION, FROM_STAR_FILE, FROM_SCIPION, \
    MEMBRANE_OUTER_SURROUNDINGS, MEMBRANE
from pyseg.protocols import *
from tomo.objects import SetOfTomoMasks, TomoMask
from tomo.protocols import ProtImportTomograms, ProtImportTomomasks


class TestFromPresegToPicking(BaseTest):
    """ """
    outputPath = None
    ds = None
    boxSize = 44
    samplingRate = 13.68

    inTomoSet = None
    inTomoSetBinned = None
    inTomomaskSetBinned = None
    protPreseg = None
    ProtGraphs = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        ds = DataSet.getDataSet('emd_10439')
        cls.ds = ds
        cls.inTomoSet = cls._importTomograms()
        cls.inTomoSetBinned = cls._normalizeTomo()
        cls.inTomomaskSetBinned = cls._ImportTomoMasks()
        cls.protPreseg = cls._runPreseg()
        cls.ProtGraphs = cls._runGraphs()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile('tomoEmd10439'),
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, 'outputTomograms', None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated.')

        return outputTomos

    @classmethod
    def _normalizeTomo(cls):
        print(magentaStr("\n==> Tomogram normalization:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=cls.inTomoSet,
                                                binning=2)

        cls.launchProtocol(protTomoNormalization)
        outputTomos = getattr(protTomoNormalization, 'outputSetOfTomograms', None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated in tomo normalization.')

        return outputTomos

    @classmethod
    def _ImportTomoMasks(cls):
        print(magentaStr("\n==> Importing data - tomoMasks"
                         ":"))
        protImportTomomasks = cls.newProtocol(ProtImportTomomasks,
                                              filesPath=cls.ds.getFile('tomomaskAnnotated'),
                                              inputTomos=cls.inTomoSetBinned)

        cls.launchProtocol(protImportTomomasks)
        tomoMaskSet = getattr(protImportTomomasks, 'outputTomoMasks', None)
        cls.assertIsNotNone(tomoMaskSet, 'No tomograms were genetated.')

        return tomoMaskSet

    @classmethod
    def _runPreseg(cls):
        print(magentaStr("\n==> Running preSeg:"))
        protPreseg = cls.newProtocol(
            ProtPySegPreSegParticles,
            segmentationFrom=FROM_SCIPION,
            inTomoMasks=cls.inTomomaskSetBinned,
            spOffVoxels=22,
            sgMembThk=60,
            sgMembNeigh=330
        )
        protPreseg.setObjLabel('Preseg')
        protPreseg = cls.launchProtocol(protPreseg)

        return protPreseg

    def testPreseg(self):
        # Check that resulting files are created as expected
        protPreseg = self.protPreseg
        nVesicles = 3
        outputStar = 'presegVesiclesCentered_pre.star'
        outputVesiclesPattern = 'emd_10439_tid_%i.mrc'
        outputVesiclesSegPattern = 'emd_10439_tid_%i_seg.mrc'
        self.assertTrue(exists(protPreseg._getExtraPath(outputStar)))
        for i in range(nVesicles):
            self.assertTrue(exists(protPreseg._getExtraPath('segs', outputVesiclesPattern % i)))
            self.assertTrue(exists(protPreseg._getExtraPath('segs', outputVesiclesSegPattern % i)))

        # Check the generated outputs
        nVesicles = 3
        binnedSamplingRate = 2 * self.samplingRate
        vesicle0Size = (186, 178, 135)
        vesicle1Size = (282, 180, 143)
        vesicle2Size = (133, 133, 139)
        vesicleSizeList = [vesicle0Size, vesicle1Size, vesicle2Size]
        setOfVesiclesTomomasks = getattr(protPreseg, 'outputTomoMasks', None)
        setOfVesiclesSubtomograms = getattr(protPreseg, 'outputSubTomograms', None)
        self.assertSetSize(setOfVesiclesTomomasks, nVesicles)
        self.assertSetSize(setOfVesiclesSubtomograms, nVesicles)
        self.assertEqual(setOfVesiclesTomomasks.getSamplingRate(), binnedSamplingRate)
        self.assertEqual(setOfVesiclesSubtomograms.getSamplingRate(), binnedSamplingRate)
        for vesicle, vesicleMask in zip(setOfVesiclesSubtomograms, setOfVesiclesTomomasks):
            self.assertEqual(vesicle.getSamplingRate(), binnedSamplingRate)
            self.assertEqual(vesicleMask.getSamplingRate(), binnedSamplingRate)
            # Preseg tid assignment to the output varies from one execution to another
            self.assertTrue(vesicle.getDimensions() in vesicleSizeList)
            self.assertTrue(vesicleMask.getDimensions() in vesicleSizeList)
        return protPreseg

    @classmethod
    def _runGraphs(cls):
        print(magentaStr("\n==> Running graphs:"))
        protGraphs = cls.newProtocol(
            ProtPySegGraphs,
            presegFrom=FROM_SCIPION,
            inSegProt=cls.protPreseg,
            maxLen=330,
            numberOfThreads=2
        )

        protGraphs.setObjLabel('Graphs')
        protGraphs = cls.launchProtocol(protGraphs)
        return protGraphs

    def testGraphs(self, protGraphs):
        # Check that resulting files are created as expected
        outputStar = 'presegVesiclesCentered_pre_mb_graph.star'
        self.assertTrue(exists(protGraphs._getExtraPath(outputStar)))
        # By default, the Disperse program intermediate results directories aren't kept
        self.assertTrue(not glob(protGraphs._getExtraPath('disperse_*')))
        return protGraphs

    # def _runFils(self, graphsProt):
    #     print(magentaStr("\n==> Running fils:"))
    #     protFils = self.newProtocol(
    #         ProtPySegFils,
    #         graphsFrom=FROM_SCIPION,
    #         inGraphsProt=graphsProt,
    #         segLabelS=MEMBRANE_OUTER_SURROUNDINGS,
    #         segLabelT=MEMBRANE,
    #         gRgEud='1 30',
    #         gRgLen='1 60',
    #         gRgSin='0 2'
    #     )
    #
    #     protFils.setObjLabel('Fils')
    #     protFils = self.launchProtocol(protFils)
    #
    #     # Check that resulting files are created as expected
    #     xmlFiles = ['mb_sources.xml', 'no_mb_targets.xml']
    #     nVesicles = 3
    #     filssFilePattern = 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_%i_fil_mb_sources_to_no_mb_targets_net'
    #     graphsFilesPerVesicle = [filssFilePattern + '.pkl',
    #                              filssFilePattern + '.vtp',
    #                              filssFilePattern + '_edges.vtp',
    #                              filssFilePattern + '_edges2.vtp',
    #                              filssFilePattern + '_skel.vtp']
    #     outputStar = 'fil_mb_sources_to_no_mb_targets_net.star'
    #
    #     [self.assertTrue(exists(protFils._getExtraPath(file))) for file in xmlFiles]
    #     self.assertTrue(exists(protFils._getExtraPath(outputStar)))
    #     for i in range(nVesicles):
    #         for file in graphsFilesPerVesicle:
    #             self.assertTrue(exists(protFils._getExtraPath(file % i)))
    #
    #     return protFils
    #
    # def _runPicking(self, filsProt, inTomoSet):
    #     print(magentaStr("\n==> Running picking:"))
    #     bSize = 30
    #     protPicking = self.newProtocol(
    #         ProtPySegPicking,
    #         inFilsProt=filsProt,
    #         inTomoSet=inTomoSet,
    #         filsFrom=FROM_SCIPION,
    #         side=MEMBRANE,
    #         boxSize=bSize
    #     )
    #
    #     protPicking.setObjLabel('Picking')
    #     protPicking = self.launchProtocol(protPicking)
    #     outputCoordinates = getattr(protPicking, 'outputCoordinates', None)
    #
    #     # Check that resulting files are created as expected
    #     xmlFile = 'mb_ext.xml'
    #     nVesicles = 3
    #     pickingFilePattern = \
    #         'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_%i_fil_mb_sources_to_no_mb_targets_net_mb_ext'
    #     graphsFilesPerVesicle = [pickingFilePattern + '.vtp',
    #                              pickingFilePattern + '_peak.vtp',
    #                              pickingFilePattern + '_surf.vtp']
    #     outputStar = 'fil_mb_sources_to_no_mb_targets_net_parts.star'
    #
    #     self.assertTrue(exists(protPicking._getExtraPath(xmlFile)))
    #     self.assertTrue(exists(protPicking._getExtraPath(outputStar)))
    #     for i in range(nVesicles):
    #         for file in graphsFilesPerVesicle:
    #             self.assertTrue(exists(protPicking._getExtraPath(file % i)))
    #
    #     # Output coordinates must have an attribute named _groupId, which stores the corresponding
    #     # vesicles index. In this case there are 3 vesicles, so it should be contained in range(2)
    #     testVesicleInd = [0, 1, 2]
    #     self.assertEqual(outputCoordinates.getBoxSize(), bSize)
    #     for coord3d in outputCoordinates:
    #         self.assertTrue(int(coord3d.getGroupId()) in testVesicleInd)
    #
    #     return protPicking
    #
    # def test_workflow(self):
    #     protPreseg = self._runPreseg()
    #     protGraphs = self._runGraphs(protPreseg)
    #     protFils = self._runFils(protGraphs)
    #     self._runPicking(protFils, protPreseg.outputSetofTomograms)
    #     # Remove generated star file
    #     if exists(self.preSegStar):
    #         remove(self.preSegStar)
