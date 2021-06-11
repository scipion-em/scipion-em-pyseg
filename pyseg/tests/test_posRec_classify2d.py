from collections import Counter
from dynamo.protocols import DynamoExtraction
from pwem.protocols import ProtImportMask
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportCoordinates3D, ProtImportTomograms
from pyseg.protocols import *
from pyseg.protocols.protocol_2d_classification import AFFINITY_PROP, CC_WITHIN_MASK, AGGLOMERATIVE, KMEANS


class TestPostRecAndClassify2d(BaseTest):
    """ """
    ds = DataSet.getDataSet('pyseg')
    samplingRate = 6.87
    boxSize = 88
    nSubtomos = 50
    minNumberOfParticles = 10

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

    def _runPosRec(self, protImportMask, protExtractSubtomo):
        print(magentaStr("\n==> Pos rec:"))
        # Pos rec
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportMask, 'outputMask', None),
        )

        protPosRec = self.launchProtocol(protPosRec)
        subtomoSet = getattr(protPosRec, 'outputSetOfSubtomogram', None)
        self.assertSetSize(subtomoSet, size=self.nSubtomos)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        return protPosRec

    def _runCl2dAffinityPropagation(self, protImportMask, protPosrec):
        print(magentaStr("\n==> 2D classification with Affinity Propagation:"))
        outClasses = 2
        protCl2d = self.newProtocol(
            ProtPySegPlaneAlignClassification,
            inputSubtomos=getattr(protPosrec, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportMask, 'outputMask', None),
            clusteringAlg=AFFINITY_PROP,
            filterSize=3,
            ccMetric=CC_WITHIN_MASK,
            apPartSizeFilter=self.minNumberOfParticles,
            apCCRefFilter=0.4
        )
        self._runClProtAndCheckResults(protPosrec, protCl2d, outClasses, clustersPostProcessed=True)

    def _runClassify2dAgglomerative(self, protImportMask, protPosrec):
        print(magentaStr("\n==> 2D classification with Agglomerative Clustering:"))
        nClasses = 3
        protCl2d = self._genCl2dProtAggOrKmeans(protImportMask, protPosrec, AGGLOMERATIVE, nClasses)
        self._runClProtAndCheckResults(protPosrec, protCl2d, nClasses)

    def _runClassify2dKmeans(self, protImportMask, protPosrec):
        print(magentaStr("\n==> 2D classification with K-means:"))
        nClasses = 4
        protCl2d = self._genCl2dProtAggOrKmeans(protImportMask, protPosrec, KMEANS, nClasses)
        self._runClProtAndCheckResults(protPosrec, protCl2d, nClasses)

    # Agglomerative and k-means clustering algorithms share parameters
    def _genCl2dProtAggOrKmeans(self, protImportMask, protPosrec, clusteringAlg, nClasses):
        return self.newProtocol(
            ProtPySegPlaneAlignClassification,
            inputSubtomos=getattr(protPosrec, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportMask, 'outputMask', None),
            clusteringAlg=clusteringAlg,
            filterSize=3,
            aggNClusters=nClasses,
        )

    def _runClProtAndCheckResults(self, protPosrec, protCl2d, nOutputClasses, clustersPostProcessed=False):
        protCl2d = self.launchProtocol(protCl2d)
        inSubtomoSet = getattr(protPosrec, 'outputSetOfSubtomogram', None)
        outSubtomoSet = getattr(protCl2d, 'outputSetOfSubtomogram', None)
        outputClasses = getattr(protCl2d, 'outputClasses', None)

        # CHECK SUBTOMO OUTPUT SET
        # Output set size must be lower or equal than the input set because some classes can have been purged
        self.assertLessEqual(outSubtomoSet.getSize(), inSubtomoSet.getSize())
        self.assertEqual(outSubtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(outSubtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))
        # Using these dataset, the classification result will be 'nOutputClasses' output classes of at least a
        # specified number of elements each (self.minNumberOfParticles)
        classesDict = dict(Counter(self._getClassIdList(outSubtomoSet)))
        self.assertFalse(any(self._getClassIdList(inSubtomoSet)))  # All classId must be 0 before classification
        self.assertEqual(len(classesDict.keys()), nOutputClasses)
        if clustersPostProcessed:
            for i in range(nOutputClasses):
                self.assertGreaterEqual(classesDict[i], self.minNumberOfParticles)

        # CHECK CLASSES OUTPUT SET
        self.assertSetSize(outputClasses, nOutputClasses)
        if clustersPostProcessed:
            for cl in outputClasses:
                self.assertGreaterEqual(cl.getSize(), self.minNumberOfParticles)

    @staticmethod
    def _getClassIdList(setObject):
        return [obj.getClassId() for obj in setObject]

    def testPosrecAndClassify2d(self):
        # Import mask
        protImportMask = self._importMask()
        # Import tomogram
        protImportTomograms = self._importTomograms()
        # Import coordinates 3D
        protImportCoords3D = self._importCoordinates3D(protImportTomograms)
        # Extract coordinates (skip CTF 3D estimation to save time)
        protExtractSubtomo = self._extractSubtomograms(protImportTomograms, protImportCoords3D)
        # Pos rec (angle randomization)
        protPosrec = self._runPosRec(protImportMask, protExtractSubtomo)
        # Classification with affinity propagation
        self._runCl2dAffinityPropagation(protImportMask, protPosrec)
        # Classification with agglomerative clustering
        self._runClassify2dAgglomerative(protImportMask, protPosrec)
        # Classification with K-means
        self._runClassify2dKmeans(protImportMask, protPosrec)