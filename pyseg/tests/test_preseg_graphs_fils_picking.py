from os import remove
from os.path import join, abspath, exists
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from relion.convert import Table
from pyseg.constants import TOMOGRAM, VESICLE, PYSEG_LABEL, MASK, FROM_STAR_FILE, FROM_SCIPION
from pyseg.protocols import *


class TestFromPresegToPicking(BaseTest):
    """ """
    ds = DataSet.getDataSet('pyseg')
    preSegStar = join(ds.getPath(), 'preseg.star')
    samplingRate = 6.87

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _genPreSegStar(self):
        """Required because pyseg's preseg requires the absolute paths contained
        in the star file"""

        # Headers for Relion's CTS star file
        preSegTable = Table(columns=[TOMOGRAM,
                                     VESICLE,
                                     PYSEG_LABEL,
                                     MASK])

        # Write star file fot 1 tomogram with 3 vesicles
        nVesicles = 3
        tomoName = abspath(self.ds.getFile('presegTomo'))
        mask = abspath(self.ds.getFile('presegMask'))
        for i in range(nVesicles):
            preSegTable.addRow(tomoName,
                               mask,
                               i + 1,
                               mask)

        preSegTable.write(self.preSegStar)

    def _runPreseg(self):
        print(magentaStr("\n==> Running preSeg:"))
        self._genPreSegStar()
        protPreseg = self.newProtocol(
            ProtPySegPreSegParticles,
            segmentationFrom=FROM_STAR_FILE,
            inStar=self.preSegStar,
            spOffVoxels=30,
            sgVoxelSize=self.samplingRate,
            sgMembThk=40,
            sgMembNeigh=200
        )

        protPreseg.setObjLabel('Preseg')
        protPreseg = self.launchProtocol(protPreseg)

        # Check that resulting files are created as expected
        nVesicles = 3
        outputStar = 'presegVesiclesCentered_pre.star'
        outputVesiclesPattern = 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_%i.mrc'
        outputVesiclesSegPattern = 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_%i_seg.mrc'
        self.assertTrue(exists(protPreseg._getExtraPath(outputStar)))
        for i in range(nVesicles):
            self.assertTrue(exists(protPreseg._getExtraPath('segs', outputVesiclesPattern % i)))
            self.assertTrue(exists(protPreseg._getExtraPath('segs', outputVesiclesSegPattern % i)))

        # Check the generated outputs
        nTomograms = 1
        nVesicles = 3
        vesiclesSize = (470, 454, 256)
        tomogramSize = (1024, 1440, 300)
        setOfVesiclesTomomasks = getattr(protPreseg, 'outputSetofTomoMasks', None)
        setOfVesiclesSubtomograms = getattr(protPreseg, 'outputSetofSubTomograms', None)
        setOfTomograms = getattr(protPreseg, 'outputSetofTomograms', None)
        self.assertSetSize(setOfVesiclesTomomasks, nVesicles)
        self.assertSetSize(setOfVesiclesSubtomograms, nVesicles)
        self.assertSetSize(setOfTomograms, nTomograms)
        self.assertEqual(setOfVesiclesTomomasks.getSamplingRate(), self.samplingRate)
        self.assertEqual(setOfVesiclesSubtomograms.getSamplingRate(), self.samplingRate)
        self.assertEqual(setOfTomograms.getSamplingRate(), self.samplingRate)
        self.assertEqual(setOfVesiclesTomomasks.getDimensions(), vesiclesSize)
        self.assertEqual(setOfVesiclesSubtomograms.getDimensions(), vesiclesSize)
        self.assertEqual(setOfTomograms.getDimensions(), tomogramSize)

        return protPreseg

    def _runGraphs(self, presegProt):
        print(magentaStr("\n==> Running graphs:"))
        protGraphs = self.newProtocol(
            ProtPySegGraphs,
            presegFrom=FROM_SCIPION,
            inSegProt=presegProt,
        )

        protGraphs.setObjLabel('Graphs')
        protGraphs = self.launchProtocol(protGraphs)

        # Check that resulting files are created as expected
        # TODO: check with Antonio which of the generated files (a lot) should be kept and add them to the test checkings
        outputStar = 'presegVesiclesCentered_pre_mb_graph.star'
        self.assertTrue(exists(protGraphs._getExtraPath(outputStar)))

        return protGraphs

    def _runFils(self, graphsProt):
        print(magentaStr("\n==> Running fils:"))
        protFils = self.newProtocol(
            ProtPySegFils,
            graphsFrom=FROM_SCIPION,
            inGraphsProt=graphsProt,
            segLabelS=3,  # From out of the membrane (labelled as 3 for this data)
            segLabelT=2,  # To inside the vesicle (labelled as 2 for this data)
            gRgEud='1 30',
            gRgLen='1 60',
            gRgSin='0 2'
        )

        protFils.setObjLabel('Fils')
        protFils = self.launchProtocol(protFils)

        # Check that resulting files are created as expected
        xmlFiles = ['mb_sources.xml', 'no_mb_targets.xml']
        nVesicles = 3
        filssFilePattern = 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_%i_fil_mb_sources_to_no_mb_targets_net'
        graphsFilesPerVesicle = [filssFilePattern + '.pkl',
                                 filssFilePattern + '.vtp',
                                 filssFilePattern + '_edges.vtp',
                                 filssFilePattern + '_edges2.vtp',
                                 filssFilePattern + '_skel.vtp']
        outputStar = 'fil_mb_sources_to_no_mb_targets_net.star'

        [self.assertTrue(exists(protFils._getExtraPath(file))) for file in xmlFiles]
        self.assertTrue(exists(protFils._getExtraPath(outputStar)))
        for i in range(nVesicles):
            for file in graphsFilesPerVesicle:
                self.assertTrue(exists(protFils._getExtraPath(file % i)))

        return protFils

    def _runPicking(self, filsProt, inTomoSet):
        print(magentaStr("\n==> Running picking:"))
        bSize = 30
        protPicking = self.newProtocol(
            ProtPySegPicking,
            inFilsProt=filsProt,
            inTomoSet=inTomoSet,
            filsFrom=FROM_SCIPION,
            side=3,  # Pick out of the membrane, labelled as 3 for this data
            boxSize=bSize
        )

        protPicking.setObjLabel('Picking')
        protPicking = self.launchProtocol(protPicking)
        outputCoordinates = getattr(protPicking, 'outputCoordinates', None)
        outputTomograms = getattr(protPicking, 'outputTomograms', None)

        # Check that resulting files are created as expected
        xmlFile = 'mb_ext.xml'
        nVesicles = 3
        pickingFilePattern = \
            'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_%i_fil_mb_sources_to_no_mb_targets_net_mb_ext'
        graphsFilesPerVesicle = [pickingFilePattern + '.vtp',
                                 pickingFilePattern + '_peak.vtp',
                                 pickingFilePattern + '_surf.vtp']
        outputStar = 'fil_mb_sources_to_no_mb_targets_net_parts.star'

        self.assertTrue(exists(protPicking._getExtraPath(xmlFile)))
        self.assertTrue(exists(protPicking._getExtraPath(outputStar)))
        for i in range(nVesicles):
            for file in graphsFilesPerVesicle:
                self.assertTrue(exists(protPicking._getExtraPath(file % i)))

        # Output coordinates must have an attribute named _pysegMembrane, which stores the corresponding
        # membrane file
        self.assertEqual(outputCoordinates.getBoxSize(), bSize)
        for coord3d in outputCoordinates:
            self.assertTrue(exists(coord3d._pysegMembrane.get()))
        # Validate precedents
        self.assertEqual(outputCoordinates.getPrecedents(), outputTomograms)

        # Validate output tomograms
        self.assertSetSize(outputTomograms, size=1)
        self.assertEqual(outputTomograms.getSamplingRate(), self.samplingRate)
        self.assertEqual(outputTomograms.getDim(), (1024, 1440, 300))

        return protPicking

    def test_workflow(self):
        protPreseg = self._runPreseg()
        protGraphs = self._runGraphs(protPreseg)
        protFils = self._runFils(protGraphs)
        self._runPicking(protFils, protPreseg.outputSetofTomograms)
        # Remove generated star file
        if exists(self.preSegStar):
            remove(self.preSegStar)
