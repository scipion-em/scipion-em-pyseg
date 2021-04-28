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
from os.path import abspath

from pwem.protocols import EMProtocol, PointerParam
from pyworkflow import BETA
from pyworkflow.protocol import EnumParam, IntParam, LEVEL_ADVANCED, GT, FloatParam, GE, LT, BooleanParam
from pyworkflow.utils import Message, makePath
from reliontomo.convert import writeSetOfSubtomograms
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import PLANE_ALIGN_CLASS_OUT, PLANE_ALIGN_CLASS_SCRIPT

# Processing level choices
PARTICLE_FLATENNING = 0     # Particle flattening
CC_MATRIX_FEAT_VECTORS = 1  # Cross correlation matrix / Feature vectors
FULL_CLASSIFICATION = 2     # Full classification

# Distance metric choices
CC_DISTANCE = 0             # Cross correlation used as distance
EUCLIDEAN_DISTANCE = 1      # Euclidean distance among image vectors used as distance metric.

# Cross correlation metric choices
CC_WITHIN_MASK = 0          # Cross-correlation within mask
SIMILARITY = 1              # Mask normalized similarity (negative squared Euclidean distance)
FULL_CC = 2                 # Full cross-correlation

# Clustering algorithm choices
AFFINITY_PROP = 0           # Affinity propagation
AGGLOMERATIVE = 1           # Agglomerative clustering
KMEANS = 2                  # K-means

# Affinity propagation condition
AP_CONDITION = 'clusteringAlg == %s and procLevel == %s' % (AFFINITY_PROP, FULL_CLASSIFICATION)

# Affinity Propagation reference choices
EXEMPLAR = 0
AVERAGE = 1


class ProtPySegPlaneAlignClassification(EMProtocol, ProtTomoBase):
    """Unsupervised and deterministic classification of membrane-bound particles"""

    _label = '2D classification'
    _devStatus = BETA
    inStarName = 'input_particles.star'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label="Input subtomograms",
                      help='Select the input subtomograms desired to be classified.')
        form.addParam('inMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask',
                      important=True,
                      allowsNull=False,
                      help='Mask used for the post processing')
        form.addParam('filterSize', IntParam,
                      label='Filter size (voxels)',
                      important=True,
                      allowsNull=False,
                      help='A value of n means that the voxels will be grouped in groups of size n.'
                      )
        form.addParam('procLevel', EnumParam,
                      label='Processing level',
                      choices=['Particle flattening',
                               'Cross correlation Z-axis radial averages',
                               'Full classification'],
                      default=FULL_CLASSIFICATION)

        group = form.addGroup('Pre-processing',)
        group.addParam('doCC3d', BooleanParam,
                       label='Do 3D radial compensation?',
                       default=True,
                       help='If "No" is selected, the normalized corss correlation (NCC) is made in 2D, '
                            'otherwise radial average is compensated for doing a NCC in 3D.')

        group = form.addGroup('Radial averages', condition='procLevel >= %s' % CC_MATRIX_FEAT_VECTORS)
        group.addParam('distanceMetric', EnumParam,
                       label='Distance metric',
                       choices=['Cross correlation',
                                'Euclidean distance among image vectors'],
                       default=CC_DISTANCE)
        group.addParam('ccMetric', EnumParam,
                       label='Cross correlation metric',
                       choices=['Cross-correlation within the mask',
                                'Mask normalized similarity',
                                'Full cross-correlation'],
                       default=CC_WITHIN_MASK,
                       condition='distanceMetric == %s' % CC_DISTANCE,
                       help='Metric used when computing the cross correlation matrix among the 2D particles. '
                            'Considerations:\n'
                            '\t- Mask normalized similarity is referred to negative squared Euclidean distance.\n'
                            '\t- Full cross-correlation (slower than ross-correlation within the mask but allows '
                            'small disalignments between particles).')
        group.addParam('pcaComps', IntParam,
                       label='PCA components for dim. reduction',
                       validators=[GT(0)],
                       default=20,
                       condition='distanceMetric == %s' % EUCLIDEAN_DISTANCE,
                       help='Number of components (moments) after the reductions.\nIf None, '
                            'then n_comp == min(n_samples, n_features) - 1')

        group = form.addGroup('Classification', condition='procLevel >= %s' % FULL_CLASSIFICATION)
        group.addParam('clusteringAlg', EnumParam,
                       label='Clustering algorithm',
                       choices=['Affinity propagation','Agglomerative clustering', 'K-means'],
                       default=AFFINITY_PROP)
        group.addParam('aggNClusters', IntParam,
                       label='Number of clusters to find',
                       default=50,
                       condition='clusteringAlg in [%s, %s]' % (AGGLOMERATIVE, KMEANS))

        group = form.addGroup('Affinity propagation',
                              condition=AP_CONDITION)
        group.addParam('apPref', FloatParam,
                       label='Affinity propagation preference (-inf, inf)',
                       default=-6,
                       help='Preference parameter (-inf, inf).\nThe smaller value the higher number of '
                            'potential classes.\nIf None, the median of the affinity class is considered.')
        group.addParam('apDumping', FloatParam,
                       label='Dumping [0.5, 1)',
                       default=0.5,
                       validators=[GE(0.5), LT(1)],
                       expertLevel=LEVEL_ADVANCED,
                       help='Dumping parameter [0.5, 1), it controls convergence speed.')
        group.addParam('apMaxIter', IntParam,
                       label='Maximum number of iterations',
                       default=2000,
                       expertLevel=LEVEL_ADVANCED)
        group.addParam('apConvIter', IntParam,
                       label='Iterations for fitting the convergence criteria',
                       default=40,
                       expertLevel=LEVEL_ADVANCED)
        group.addParam('apReference', EnumParam,
                       label='Reference 2D image used for classes',
                       choices=['Exemplar', 'Average'],
                       default=EXEMPLAR,
                       expertLevel=LEVEL_ADVANCED)

        group = form.addGroup('Classification post-processing',
                              condition=AP_CONDITION,
                              expertLevel=LEVEL_ADVANCED)
        group.addParam('apPartSizeFilter', IntParam,
                       label='Minimum number of particles per class',
                       validators=[GE(0)],
                       default=0,
                       expertLevel=LEVEL_ADVANCED,
                       help='Purge classes with less than the specified number of particles.')
        group.addParam('apCCRefFilter', FloatParam,
                       label='Cross-correlation against AP reference filter',
                       validators=[GE(0)],
                       default=0,
                       expertLevel=LEVEL_ADVANCED,
                       help='Purge classes with the cross correlation against the reference '
                            'lower than the specified value.')
        form.addParallelSection(threads=4, mpi=0)
        # group.addParam('aggLinkage', EnumParam,
        #                label='linkage criterion',
        #                choices=[])
        # TODO: ask Antonio about the linkage criterion, if other values apart from ward are implemented (it seems not
        # to be the case)

    def _insertAllSteps(self):
        # self._initialize()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('pysegPlaneAlignClassification')
        # self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        subtomoSet = self.inputSubtomos.get()
        subTomoStar = self._getExtraPath(self.inStarName)
        writeSetOfSubtomograms(subtomoSet, subTomoStar, isPyseg=True)

    def pysegPlaneAlignClassification(self):
        # Generate output subtomo dir
        outDir = self._getExtraPath(PLANE_ALIGN_CLASS_OUT)
        makePath(outDir)
        # Script called
        Plugin.runPySeg(self, PYTHON, self. _getCommand(outDir))

    # def createOutputStep(self, outStar):
    #     # Read generated star file and create the output objects
    #     self.subtomoSet = self._createSetOfSubTomograms()
    #     self.subtomoSet.setSamplingRate(self.inMask.get().getSamplingRate())
    #     warningMsg = readStarFile(self, self.subtomoSet, RELION_SUBTOMO_STAR, starFile=outStar)
    #     if warningMsg:
    #         self.warningMsg = String(warningMsg)
    #         self._store()
    #
    #     self._defineOutputs(outputSetOfSubtomogram=self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            procLevel = self.procLevel.get()
            if procLevel == PARTICLE_FLATENNING:
                summary.append('Processing level: *Particles flattening*')
            elif procLevel == CC_MATRIX_FEAT_VECTORS:
                summary.append('Processing level: *CC matrix or Feature vectors computing*')
            if procLevel >= PARTICLE_FLATENNING:
                r3dMsg = 'Radial compensation for 3D'
                r3dMsg = r3dMsg if self.doCC3d.get() else 'No ' + r3dMsg
                summary.append(
                    '\nParticles pre-processing:\n'
                    '   - Low pass Gaussian filter sigma: %i voxels.\n'
                    '   - %s' % (self.filterSize.get(), r3dMsg)
                )
            if procLevel >= CC_MATRIX_FEAT_VECTORS:
                distanceMetric = self.distanceMetric.get()
                if distanceMetric == CC_DISTANCE:
                    msg = '   - Cross correlation metric: %s\n' % self._decodeCCMetric()
                else:
                    msg = '   - Number of components for PCA dimensionality reduction: %i\n' % self.pcaComps.get()
                summary.append(
                    '\nDistance metric calculation:\n'
                    '   - Distance metric: %s\n%s' % (self._decodeDistanceMetric(), msg))
        return summary

    def _validate(self):
        alg = self.clusteringAlg.get()
        distMetric = self.distanceMetric.get()
        nSubtomos = self.inputSubtomos.get().getSize()
        errors = []
        if distMetric == EUCLIDEAN_DISTANCE and self.pcaComps.get() > nSubtomos:
            errors.append('Number of PCA components must be between 0 and min(n_samples, n_features).')
        if alg == KMEANS and distMetric != EUCLIDEAN_DISTANCE:
            errors.append('*K-means* clustering algorithm requires *Euclidean distance* as distance metric.')
        if alg != AFFINITY_PROP and self.aggNClusters.get() <= 0 or self.aggNClusters.get() > nSubtomos:
            errors.append('Number of clusters to find must be in range (0, nParticles).')
        return errors

    # --------------------------- UTIL functions -----------------------------------

    # def _initialize(self):
    #     # Only the ecludian distance is a valid distance metric if the clustering algorithm is different to
    #     # Affinity Propagation
    #     # if self.clusteringAlg.get() != AFFINITY_PROP:
    #     #     self.distanceMetric.set(EUCLIDEAN_DISTANCE)
    #     if self.distanceMetric.get() == CC_DISTANCE:
    #         self.ccMetric.set()


    def _getCommand(self, outDir):
        alg = self.clusteringAlg.get()
        classCmd = ' '
        classCmd += '%s ' % Plugin.getHome(PLANE_ALIGN_CLASS_SCRIPT)
        classCmd += '--inRootDir scipion '
        classCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        classCmd += '--inMask %s ' % abspath(self.inMask.get().getFileName())
        classCmd += '--outDir %s ' % outDir
        classCmd += '--filterSize %s ' % self.filterSize.get()
        classCmd += '--procLevel %s ' % (self.procLevel.get() + 1)  # Numbered from 1 in pyseg
        classCmd += '--doCC3d %s ' % self.doCC3d.get()
        classCmd += '--ccMetric %s ' % self._decodeCCMetric()
        classCmd += '--clusteringAlg %s ' % self._decodeClusteringAlg()
        classCmd += '--distanceMetric %s ' % self._decodeDistanceMetric()
        classCmd += '--pcaComps %s ' % self.pcaComps.get()
        if alg == AFFINITY_PROP:
            classCmd += '--apPref %s ' % self.apPref.get()
            classCmd += '--apDumping %s ' % self.apDumping.get()
            classCmd += '--apMaxIter %s ' % self.apMaxIter.get()
            classCmd += '--apConvIter %s ' % self.apConvIter.get()
            classCmd += '--apPartSizeFilter %s ' % self.apPartSizeFilter.get()
            classCmd += '--apCCRefFilter %s ' % self.apCCRefFilter.get()
        elif alg == AGGLOMERATIVE:
            classCmd += '--aggNClusters %s ' % self.aggNClusters.get()
        else:
            classCmd += '--kmeansNClusters %s ' % self.aggNClusters.get()

        classCmd += '-j %s ' % self.numberOfThreads.get()

        return classCmd

    def _decodeCCMetric(self):
        res = None
        ccMetric = self.ccMetric.get()
        if ccMetric == CC_WITHIN_MASK:
            res = 'cc'
        elif ccMetric == SIMILARITY:
            res = 'similarity'
        elif ccMetric == FULL_CC:
            res = 'full_cc'

        return res

    def _decodeClusteringAlg(self):
        res = None
        clusteringAlg = self.clusteringAlg.get()
        if clusteringAlg == AFFINITY_PROP:
            res = 'AP'
        elif clusteringAlg == AGGLOMERATIVE:
            res = 'AG'
        elif clusteringAlg == KMEANS:
            res = 'Kmeans'

        return res

    def _decodeDistanceMetric(self):
        res = None
        distanceMetric = self.distanceMetric.get()
        if distanceMetric == CC_DISTANCE:
            res = 'ncc_2dz'
        elif distanceMetric == EUCLIDEAN_DISTANCE:
            res = 'vectors'

        return res

