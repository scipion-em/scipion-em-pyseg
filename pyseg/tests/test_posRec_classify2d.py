# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************

from collections import Counter
from xmipp3.constants import MASK3D_CYLINDER
from xmipp3.protocols import XmippProtCreateMask3D
from xmipp3.protocols.protocol_preprocess.protocol_create_mask3d import SOURCE_GEOMETRY

from pyseg.protocols import ProtPySegPostRecParticles
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols import ProtImportSubtomogramsFromStar
from pyseg.protocols.protocol_2d_classification import AFFINITY_PROP, CC_WITHIN_MASK, AGGLOMERATIVE, KMEANS, \
    ProtPySegPlaneAlignClassification
from pyseg.protocols.protocol_2d_classification import outputObjects as cl2dOutputs
from reliontomo.protocols.protocol_import_subtomograms_from_star import outputObjects as importSubtomoOutputs
from tomo.protocols import ProtImportTomograms
from tomo.tests import EMD_10439, DataSetEmd10439


class TestPostRecAndClassify2d(BaseTest):

    dataset = None
    samplingRate = 13.68
    boxSize = 44
    tomoId = 'emd_10439'
    nSubtomos = 7
    minNumberOfParticles = 2
    protImporSubtomogramsFromStar = None
    protCreateParticleMask = None
    protCreateMembraneMask = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet(EMD_10439)
        inTomoSet = cls._importTomograms()
        cls.protImporSubtomogramsFromStar = cls._runImportSubtomogramsFromStarFile(inTomos=inTomoSet)
        cls.protCreateParticleMask = cls._runCreate3dParticleMask()
        cls.protCreateMembraneMask = cls._runCreate3dMembraneMask()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.dataset.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, 'outputTomograms', None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated.')
        return outputTomos

    @classmethod
    def _runImportSubtomogramsFromStarFile(cls, inTomos=None):
        protImporSubtomogramsFromStar = cls.newProtocol(ProtImportSubtomogramsFromStar,
                                                        starFile=cls.dataset.getFile(DataSetEmd10439.subtomogramsStarFile.name),
                                                        inTomos=inTomos,
                                                        samplingRate=cls.samplingRate,
                                                        boxSize=cls.boxSize)

        return cls.launchProtocol(protImporSubtomogramsFromStar)

    @classmethod
    def _runCreate3dMask(cls, doShift=False, zShift=0, height=20, doSmooth=False):
        protMask3D = cls.newProtocol(XmippProtCreateMask3D,
                                     source=SOURCE_GEOMETRY,
                                     samplingRate=cls.samplingRate,
                                     size=cls.boxSize,
                                     geo=MASK3D_CYLINDER,
                                     radius=12,
                                     shiftCenter=doShift,
                                     centerZ=zShift,
                                     height=height,
                                     doSmooth=doSmooth
                                     )
        return cls.launchProtocol(protMask3D)

    @classmethod
    def _runCreate3dParticleMask(cls):
        return cls._runCreate3dMask(doShift=True, zShift=6, height=30, doSmooth=True)

    @classmethod
    def _runCreate3dMembraneMask(cls):
        return cls._runCreate3dMask(doShift=False, height=5)

    def testPosRecOnlyRotRand(self):
        print(magentaStr("\n==> Pos rec with rot angle randomization only:"))
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(self.protImporSubtomogramsFromStar, importSubtomoOutputs.subtomograms.name),
            inMask=getattr(self.protCreateParticleMask, 'outputMask', None),
        )
        protPosRec.setObjLabel('Pos rec only rot rand')
        protPosRec = self.launchProtocol(protPosRec)
        self._runCheckPosRec(protPosRec)
        return protPosRec

    def testPosRecWithMbAttenuation(self):
        print(magentaStr("\n==> Pos rec with membrane attenuation:"))
        protPosRec = self.newProtocol(
            ProtPySegPostRecParticles,
            inputSubtomos=getattr(self.protImporSubtomogramsFromStar, importSubtomoOutputs.subtomograms.name),
            inMask=getattr(self.protCreateParticleMask, 'outputMask', None),
            mbMask=getattr(self.protCreateMembraneMask, 'outputMask', None),
            mbSupFactor=0.3
        )
        protPosRec.setObjLabel('Pos rec with membrane attenuation')
        protPosRec = self.launchProtocol(protPosRec)
        self._runCheckPosRec(protPosRec)
        return protPosRec

    def _runCheckPosRec(self, protPosRec):
        subtomoSet = getattr(protPosRec, cl2dOutputs.subtomograms.name)
        self.assertSetSize(subtomoSet, size=self.nSubtomos)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

    def testCl2dAffinityPropagation(self):
        print(magentaStr("\n==> 2D classification with Affinity Propagation:"))
        outClasses = 1
        protCl2d = self.newProtocol(
            ProtPySegPlaneAlignClassification,
            inputSubtomos=getattr(self.protImporSubtomogramsFromStar, importSubtomoOutputs.subtomograms.name),
            inMask=getattr(self.protCreateParticleMask, 'outputMask', None),
            clusteringAlg=AFFINITY_PROP,
            filterSize=2,
            ccMetric=CC_WITHIN_MASK,
            apPartSizeFilter=self.minNumberOfParticles
        )
        protCl2d.setObjLabel('cl2d - AP')
        self._runClProtAndCheckResults(protCl2d, outClasses, clustersPostProcessed=True)

    def testClassify2dAgglomerative(self):
        print(magentaStr("\n==> 2D classification with Agglomerative Clustering:"))
        nClasses = 3
        protCl2d = self._genCl2dProtAggOrKmeans(AGGLOMERATIVE, nClasses)
        protCl2d.setObjLabel('cl2d - AG')
        self._runClProtAndCheckResults(protCl2d, nClasses)

    def testClassify2dKmeans(self):
        print(magentaStr("\n==> 2D classification with K-means:"))
        nClasses = 4
        protCl2d = self._genCl2dProtAggOrKmeans( KMEANS, nClasses)
        protCl2d.setObjLabel('cl2d - K-means')
        self._runClProtAndCheckResults(protCl2d, nClasses)

    # Agglomerative and k-means clustering algorithms share parameters
    def _genCl2dProtAggOrKmeans(self, clusteringAlg, nClasses):
        return self.newProtocol(
            ProtPySegPlaneAlignClassification,
            inputSubtomos=getattr(self.protImporSubtomogramsFromStar, importSubtomoOutputs.subtomograms.name),
            inMask=getattr(self.protCreateParticleMask, 'outputMask', None),
            clusteringAlg=clusteringAlg,
            filterSize=2,
            aggNClusters=nClasses,
            pcaComps=3
        )

    def _runClProtAndCheckResults(self, protCl2d, nOutputClasses, clustersPostProcessed=False):
        protCl2d = self.launchProtocol(protCl2d)
        inSubtomoSet = getattr(self.protImporSubtomogramsFromStar, importSubtomoOutputs.subtomograms.name)
        outSubtomoSet = getattr(protCl2d, cl2dOutputs.subtomograms.name)
        outputClasses = getattr(protCl2d, cl2dOutputs.classes.name)

        # CHECK SUBTOMO OUTPUT SET
        # Output set size must be lower or equal than the input set because some classes can have been purged
        self.assertLessEqual(outSubtomoSet.getSize(), inSubtomoSet.getSize())
        self.assertEqual(outSubtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(outSubtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))
        # Using these dataset, the classification result will be 'nOutputClasses' output classes of at least a
        # specified number of elements each (self.minNumberOfParticles)
        classesDict = dict(Counter(self._getClassIdList(outSubtomoSet)))
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
