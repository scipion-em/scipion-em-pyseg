from collections import Counter
from os.path import exists

from dynamo.protocols import DynamoExtraction
from imod.protocols import ProtImodAutomaticCtfEstimation
from pwem.protocols import ProtImportMask
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols.protocol_ctf_3d_estimation import CTF3D_PER_SUBVOLUME, ProtRelionEstimateCTF3D
from tomo.protocols import ProtImportCoordinates3D, ProtImportTomograms, ProtImportTs
from pyseg.protocols import *
from pyseg.protocols.protocol_2d_classification import AFFINITY_PROP, CC_WITHIN_MASK, AGGLOMERATIVE, KMEANS


class TestPostRecAndClassify2d(BaseTest):
    """ """
    samplingRate = 6.87  # Angstrom/voxel
    boxSize = 40  # pixel
    nSubtomos = 50
    minNumberOfParticles = 10

    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('pyseg')
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('pyseg')

    def _importCylMask(self):
        print(magentaStr("\n==> Importing the cylinder mask:"))
        protImportMask = self._importMask(self.ds.getFile('posRecMask'))
        protImportMask.setObjLabel('Import cyl mask')
        return protImportMask

    def _importMembraneMask(self):
        print(magentaStr("\n==> Importing the membrane mask:"))
        protImportMask = self._importMask(self.ds.getFile('mbMask'))
        protImportMask.setObjLabel('Import membrane mask')
        return protImportMask

    def _importMask(self, maskPath):
        protImportMask = self.newProtocol(
            ProtImportMask,
            maskPath=maskPath,
            samplingRate=self.samplingRate,
        )
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

    def _importTiltSeries(self):
        print(magentaStr("\n==> Importing data - tilt series:"))
        protImportTS = self.newProtocol(
            ProtImportTs,
            filesPath=self.ds.getPath(),
            filesPattern='Pertuzumab_1_defocus_25um_{TS}_ali.mrc',
            anglesFrom=2,  # ANGLES_FROM_TLT
            voltage=300,
            magnification=10000,
            samplingRate=self.samplingRate,
            tiltAxisAngle=85.2,
            dosePerFrame=2.1
        )
        protImportTS.setObjLabel('import tilt series')
        protImportTS = self.launchProtocol(protImportTS)
        tsSet = getattr(protImportTS, 'outputTiltSeries', None)

        # Validate output tomograms
        self.assertSetSize(tsSet, size=1)
        self.assertEqual(tsSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(tsSet.getDim(), (4092, 5760, 41))
        return protImportTS

    def _tSCtfEstimateImod(self, protImportTS):
        print(magentaStr("\n==> Estimating the tilt series ctf:"))
        protTSCtfImod = self.newProtocol(
            ProtImodAutomaticCtfEstimation,
            inputSet=getattr(protImportTS, 'outputTiltSeries', None),
        )
        protTSCtfImod.setObjLabel('estimate TS CTF')
        protTSCtfImod = self.launchProtocol(protTSCtfImod)
        ctfSeriesSet = getattr(protTSCtfImod, 'outputSetOfCTFTomoSeries', None)

        # Validate output tomograms
        self.assertSetSize(ctfSeriesSet, size=1)

        return protTSCtfImod

    def _estimateCTF3D(self, protImportCoords3D, protTSCtfImod):
        print(magentaStr("\n==> Estimating the 3D CTF per subvolume:"))
        protEstimateCTF3D = self.newProtocol(
            ProtRelionEstimateCTF3D,
            inputCoordinates=getattr(protImportCoords3D, 'outputCoordinates', None),
            inputSetCTFTomoSeries=getattr(protTSCtfImod, 'outputSetOfCTFTomoSeries', None),
            doseFilesPath=self.ds.getPath(),
            filesPattern='*ExpDose.txt',
            boxSize=self.boxSize,
            ctf3dMode=CTF3D_PER_SUBVOLUME,
        )
        protEstimateCTF3D.setObjLabel('Estimate CTF 3D')
        protEstimateCTF3D = self.launchProtocol(protEstimateCTF3D)
        coord3DSet = getattr(protEstimateCTF3D, 'outputCoordinates', None)

        # Validate output tomograms
        self.assertSetSize(coord3DSet, size=self.nSubtomos)
        self.assertEqual(coord3DSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(coord3DSet.getBoxSize(), self.boxSize)
        # Output coordinates must have an attribute named _3dcftMrcFile, which stores the
        # path of the each ctf3D file
        for coord3d in coord3DSet:
            self.assertTrue(exists(coord3d._3dcftMrcFile.get()))

        return protEstimateCTF3D

    def _extractSubtomograms(self, protImportTomo, protEstimateCTF3D):
        print(magentaStr("\n==> Extracting subtomograms:"))
        protDynamoExtract = self.newProtocol(
            DynamoExtraction,
            inputCoordinates=getattr(protEstimateCTF3D, 'outputCoordinates', None),
            tomoSource=1,  # OTHER
            boxSize=self.boxSize,
            doInvert=True,
            inputTomograms=getattr(protImportTomo, 'outputTomograms', None)
        )
        protDynamoExtract = self.launchProtocol(protDynamoExtract)

        return protDynamoExtract

    def _runPosRecOnlyRotRand(self, protImportCykMask, protExtractSubtomo):
        print(magentaStr("\n==> Pos rec with rot angle randomization only:"))
        # Pos rec
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportCykMask, 'outputMask', None),
            # mbMask=getattr(protImportMembraneMask, 'outputMask', None),
            # mbSupFactor=0.3
        )
        protPosRec.setObjLabel('Pos rec only rot rand')
        protPosRec = self.launchProtocol(protPosRec)
        self._runCheckPosRec(protPosRec)
        return protPosRec

    def _runPosRecWithMbAttenuation(self, protImportCykMask, protImportMembraneMask, protExtractSubtomo):
        print(magentaStr("\n==> Pos rec with membrane attenuation:"))
        # Pos rec
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportCykMask, 'outputMask', None),
            mbMask=getattr(protImportMembraneMask, 'outputMask', None),
            mbSupFactor=0.3
        )
        protPosRec.setObjLabel('Pos rec with membrane attenuation')
        protPosRec = self.launchProtocol(protPosRec)
        self._runCheckPosRec(protPosRec)
        return protPosRec

    def _runPosRecWithMbAttenuationAndFilter(self, protImportCykMask, protImportMembraneMask, protExtractSubtomo):
        print(magentaStr("\n==> Pos rec with membrane attenuation and applying gaussian low pass filter:"))
        # Pos rec
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            inMask=getattr(protImportCykMask, 'outputMask', None),
            mbMask=getattr(protImportMembraneMask, 'outputMask', None),
            mbSupFactor=0.3,
            doGaussLowPassFilter=True,
            cutOffRes=1,
            ampCutOff=0.01,
            filterCTF=True
        )
        protPosRec.setObjLabel('Pos rec with mb att and glpf')
        protPosRec = self.launchProtocol(protPosRec)
        self._runCheckPosRec(protPosRec)
        return protPosRec

    def _runCheckPosRec(self, protPosRec):
        subtomoSet = getattr(protPosRec, 'outputSetOfSubtomogram', None)
        self.assertSetSize(subtomoSet, size=self.nSubtomos)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

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

    def _runCommonPrevProtocolsRequired(self):
        try:
            # Import masks
            self.protImportCylMask = self._importCylMask()
            self.protImportMbMask = self._importMembraneMask()
            # Import tomogram
            self.protImportTomograms = self._importTomograms()
            # Import coordinates 3D
            self.protImportCoords3D = self._importCoordinates3D(self.protImportTomograms)
            # Import tilt series
            self.protImportTS = self._importTiltSeries()
            # Estimate TS CTF
            self.protEstimateCTF = self._tSCtfEstimateImod(self.protImportTS)
            # Estimate CTF 3D per subvolume
            self.protEstimateCTF3D = self._estimateCTF3D(self.protImportCoords3D, self.protEstimateCTF)
            # Extract coordinates
            self.protExtractSubtomo = self._extractSubtomograms(self.protImportTomograms, self.protEstimateCTF3D)
        except:
            raise 'One of the previous protocols required failed'

    def testPosRecAndClassiftWorkflow(self):
        self._runCommonPrevProtocolsRequired()
        # Posrec with only rot angle randomization
        ProtPosRecRotRand = self._runPosRecOnlyRotRand(self.protImportCylMask, self.protExtractSubtomo)
        # 2D classification with affinity propagation
        self._runCl2dAffinityPropagation(self.protImportCylMask, ProtPosRecRotRand)
        # Posrec with membrane attenuation
        ProtPosRecMbAtt = self._runPosRecWithMbAttenuation(
            self.protImportCylMask, self.protImportMbMask, self.protExtractSubtomo)
        # 2D classification with agglomerative clustering
        self._runClassify2dAgglomerative(self.protImportCylMask, ProtPosRecMbAtt)
        # Posrec with membrane attenuation and applying a gaussian low pass filter
        ProtPosRecMbAttAndFilter = self._runPosRecWithMbAttenuationAndFilter(
            self.protImportCylMask, self.protImportMbMask, self.protExtractSubtomo)
        # 2D classification with k-means
        self._runClassify2dAgglomerative(self.protImportCylMask, ProtPosRecMbAttAndFilter)