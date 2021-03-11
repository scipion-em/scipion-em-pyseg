from dynamo.protocols import DynamoExtraction
from pwem.protocols import ProtImportMask
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportCoordinates3D, ProtImportTomograms
from pyseg.protocols import *


class TestPostRec(BaseTest):
    """ """
    ds = DataSet.getDataSet('pyseg')
    samplingRate = 6.87
    boxSize = 88

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _importMask(self):
        print(magentaStr("\n==> Importing the mask:"))
        protImportMask = self.newProtocol(
            ProtImportMask,
            maskPath=self.ds.getFile('posRecMask'),
            samplingRate=self.samplingRate,
        )
        protImportMask.setObjLabel('Import mask')
        protImportMask = self.launchProtocol(protImportMask)
        outVol = getattr(protImportMask, 'outputMask', None)

        # Validate output tomograms
        self.assertEqual(outVol.getDim(), (self.boxSize, self.boxSize, self.boxSize))
        self.assertEqual(outVol.getSamplingRate(), self.samplingRate)

        return protImportMask

    def _importTomograms(self):
        print(magentaStr("\n==> Import tomograms:"))
        protImportTomograms = self.newProtocol(
            ProtImportTomograms,
            filesPath=self.ds.getFile('presegTomo'),
            samplingRate=self.samplingRate
        )
        self.launchProtocol(protImportTomograms)

        return protImportTomograms

    def _importCoordinates3D(self, protImportTomo):
        print(magentaStr("\n==> Import 3D coordinates:"))
        protImportCoords3D = self.newProtocol(
            ProtImportCoordinates3D,
            importFrom=3,  # IMPORT_FROM_DYNAMO
            filesPath=self.ds.getPath(),
            filesPattern='*.tbl',
            samplingRate=self.samplingRate,
            boxSize=self.boxSize,
            importTomograms=getattr(protImportTomo, 'outputTomograms', None)
        )
        self.launchProtocol(protImportCoords3D)

        return protImportCoords3D

    def _extractSubtomograms(self, protImportTomo, protImportCoords3D):
        print(magentaStr("\n==> Extracting subtomograms:"))
        protDynamoExtract = self.newProtocol(
            DynamoExtraction,
            inputCoordinates=getattr(protImportCoords3D, 'outputCoordinates', None),
            tomoSource=1,  # OTHER
            boxSize=self.boxSize,
            doInvert=True,
            inputTomograms=getattr(protImportTomo, 'outputTomograms', None)
        )
        protDynamoExtract = self.launchProtocol(protDynamoExtract)

        return protDynamoExtract

    def testPostRec(self):
        # Import mask
        protImportMask = self._importMask()
        # Import tomogram
        protImportTomograms = self._importTomograms()
        # Import coordinates 3D
        protImportCoords3D = self._importCoordinates3D(protImportTomograms)
        # Extract coordinates (skip CTF 3D estimation to save time)
        protExtractSubtomo = self._extractSubtomograms(protImportTomograms, protImportCoords3D)
        print(magentaStr("\n==> Pos rec:"))
        # Pos rec
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportMask, 'outputMask', None),
        )

        protPosRec = self.launchProtocol(protPosRec)

        # Check results
        nSubtomos = 6
        subtomoSet = getattr(protPosRec, 'outputSetOfSubtomogram', None)
        self.assertSetSize(subtomoSet, size=nSubtomos)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))
