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
from enum import Enum
from os.path import abspath, join, exists
import glob
from emtable import Table
from pwem.protocols import EMProtocol, PointerParam
from pyseg.convert import readPysegSubtomograms
from pyseg.utils import checkMaskFormat, getFinalMaskFileName
from pyworkflow.object import String
from pyworkflow.protocol import EnumParam, IntParam, LEVEL_ADVANCED, FloatParam, GE, LT, BooleanParam
from pyworkflow.utils import Message, makePath
from reliontomo.convert import createWriterTomo
from scipion.constants import PYTHON
from tomo.objects import SetOfSubTomograms, SetOfClassesSubTomograms
from tomo.protocols import ProtTomoBase
from pyseg import Plugin
from pyseg.constants import PLANE_ALIGN_CLASS_OUT, PLANE_ALIGN_CLASS_SCRIPT, SEE_METHODS_TAB


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
AP_CONDITION = 'clusteringAlg == %s' % AFFINITY_PROP
NOT_AP_CONDITION = 'clusteringAlg != %s' % AFFINITY_PROP

# Affinity Propagation reference choices
EXEMPLAR = 0
AVERAGE = 1


class outputObjects(Enum):
    subtomograms = SetOfSubTomograms
    classes = SetOfClassesSubTomograms


class ProtPySegPlaneAlignClassification(EMProtocol, ProtTomoBase):
    """
    Performs unsupervised classification of membrane-bound subtomograms
    by comparing their structural organization through plane-aligned
    particle representations. The protocol is designed to identify
    biologically meaningful particle groups in cryo-electron tomography
    datasets while reducing user bias and preserving structural diversity.

    AI Generated:

    Plane Align Classification (ProtPySegPlaneAlignClassification) - User Manual
        Overview

        The Plane Align Classification protocol performs unsupervised
        classification of membrane-associated subtomograms using
        rotationally aligned particle representations and similarity-based
        clustering strategies. Its main purpose is to separate heterogeneous
        particle populations into structurally related classes that can be
        interpreted biologically or used in downstream subtomogram averaging
        workflows.

        In cryo-electron tomography studies, membrane-bound complexes often
        display substantial variability in orientation, conformation,
        composition, or local membrane context. This protocol addresses
        these challenges by generating normalized particle representations
        that emphasize structural similarities while reducing irrelevant
        variability. The resulting classes can reveal distinct molecular
        states, assembly intermediates, or spatially organized membrane
        populations.

        Inputs and Biological Context

        The protocol requires a set of subtomograms together with a 3D mask.
        The subtomograms correspond to extracted particles from tomographic
        volumes, while the mask defines the region of interest that should
        contribute to the classification process. In most biological
        applications, the mask should encompass the structurally relevant
        region while excluding unnecessary solvent or noisy membrane areas.

        Careful preparation of the input data is important for meaningful
        classification. The subtomograms and mask should share the same box
        dimensions and voxel size so that particle features are represented
        consistently. Significant mismatches between masks and particles may
        lead to unstable classifications or biologically misleading groups.

        In practical cryo-ET workflows, this protocol is commonly applied
        after particle extraction and membrane orientation assignment. It is
        especially useful for studying membrane protein assemblies,
        vesicular transport systems, cytoskeletal interactions, or crowded
        membrane environments where structural heterogeneity is expected.

        Particle Pre-processing

        Before classification, the protocol performs several preprocessing
        operations intended to improve the robustness of similarity
        estimation. Filtering helps suppress high-frequency noise while
        preserving the dominant structural features of the particles. This
        is particularly important in cryo-electron tomography datasets,
        where the signal-to-noise ratio is often low.

        The filter size determines the level of smoothing applied to the
        particle representations. Small filter values preserve fine details
        but may remain sensitive to noise, whereas larger values emphasize
        global structural organization at the expense of local features.
        Biological users should adapt this parameter according to the size
        and expected flexibility of the target complex.

        The protocol also supports optional radial compensation during
        similarity estimation. This operation attempts to compensate for
        intensity biases associated with radial averaging and can improve
        comparisons between particles that exhibit strong three-dimensional
        structural variation. In many membrane-associated systems, enabling
        this option improves classification stability, especially for
        complexes extending away from the membrane surface.

        Similarity Metrics and Structural Comparison

        A key aspect of the protocol is the calculation of similarity
        relationships between particles. Different cross-correlation metrics
        are available depending on the biological problem and the desired
        balance between robustness and computational cost.

        Cross-correlation within the mask focuses the comparison on the
        masked structural region and is generally suitable for most
        membrane-associated particles. This approach minimizes the influence
        of solvent and unrelated densities. Full cross-correlation provides
        greater tolerance to small particle misalignments but may increase
        computational complexity. The similarity metric based on normalized
        distances can be useful when emphasizing overall structural
        resemblance rather than precise density overlap.

        From a biological perspective, the selected metric influences which
        structural features dominate the classification. Metrics focused on
        local density agreement tend to separate subtle conformational
        differences, whereas broader similarity measures may emphasize
        large-scale organizational patterns.

        Clustering Algorithms

        The protocol provides multiple clustering strategies suitable for
        different levels of heterogeneity and dataset complexity.

        Affinity propagation automatically determines representative classes
        from the particle similarity relationships. This method is
        particularly attractive for exploratory analyses because it does not
        require the user to define the number of classes in advance.
        Biologically, it is useful when the diversity of conformational or
        compositional states is unknown.

        Agglomerative clustering progressively groups particles according to
        structural similarity. This approach is often suitable when the user
        expects hierarchical relationships between classes or wishes to
        impose a predefined number of groups.

        K-means clustering partitions particles into a fixed number of
        classes defined by the user. It is computationally efficient and can
        perform well when the approximate heterogeneity level is already
        known from prior biological knowledge or exploratory analyses.

        Dimensionality Reduction

        For clustering methods based on feature vectors, the protocol can
        perform dimensionality reduction using principal component analysis.
        This operation compresses the structural information into a reduced
        set of components while preserving the dominant variability present
        in the dataset.

        Dimensionality reduction becomes particularly important for large
        subtomogram datasets because it improves computational efficiency
        and reduces sensitivity to noise. In biological applications, using
        too few components may oversimplify structural variability, whereas
        too many components may reintroduce noise and unstable features.

        Post-processing and Class Filtering

        The protocol includes optional post-processing operations that help
        improve the interpretability of the final classification results.
        Small classes containing very few particles can be removed to reduce
        the impact of poorly supported structural groups. This is often
        useful when analyzing noisy or highly heterogeneous cryo-ET data.

        Additional filtering based on similarity against representative
        class references can also be applied. This step helps eliminate weak
        or poorly defined classes that may arise from unstable clustering.
        Biologically, these filters should be used carefully because overly
        aggressive thresholds may discard rare but meaningful particle
        states.

        Outputs and Biological Interpretation

        After execution, the protocol produces a classified set of
        subtomograms together with the corresponding classes and
        representative reference images. Each particle is assigned to a
        structural class, allowing users to inspect particle organization,
        compare class populations, and identify biologically meaningful
        patterns.

        The representative images associated with each class summarize the
        dominant structural features of that particle population. These
        representatives can be used for visual inspection, downstream
        subtomogram averaging, or further structural refinement.

        In biological studies, the resulting classes may correspond to
        distinct conformational states, oligomeric assemblies, interaction
        partners, or membrane-associated functional states. However,
        interpretation should always consider the possibility of continuous
        heterogeneity, missing-wedge artifacts, and limited particle counts.

        Practical Recommendations

        For exploratory analyses, affinity propagation is often a good
        starting point because it adapts naturally to unknown heterogeneity.
        When the approximate number of structural states is known,
        agglomerative clustering or K-means may provide more controlled and
        reproducible partitions.

        Appropriate masking is one of the most important factors for
        successful classification. The mask should isolate the biologically
        relevant region while minimizing unrelated membrane or solvent
        signal. Poor masking frequently produces unstable or biologically
        ambiguous classes.

        Moderate filtering is usually beneficial for noisy cryo-ET data, but
        excessive smoothing may obscure subtle conformational differences.
        Users are encouraged to visually inspect representative classes and
        compare multiple parameter combinations when studying highly dynamic
        systems.

        Final Perspective

        Classification of membrane-bound subtomograms is not merely a
        computational grouping procedure but a biologically meaningful
        strategy for exploring molecular diversity inside native cellular
        environments. Reliable results depend on thoughtful preprocessing,
        biologically informed masking, and careful interpretation of class
        heterogeneity. When applied appropriately, this protocol provides a
        powerful framework for discovering structural organization and
        functional variability in cryo-electron tomography datasets.
    """

    _label = '2D classification'
    inStarName = 'input_particles.star'
    _dataTable = Table()
    _outDir = None
    _warningMsg = String()
    _distanceMetric = None

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
        form.addParam('clusteringAlg', EnumParam,
                      label='Clustering algorithm',
                      choices=['Affinity propagation', 'Agglomerative clustering', 'K-means'],
                      important=True,
                      default=AFFINITY_PROP)

        group = form.addGroup('Pre-processing',)
        group.addParam('filterSize', IntParam,
                       label='Filter size (voxels)',
                       allowsNull=False,
                       important=True,
                       help='A value of n means that the voxels will be grouped in groups of size n.'
                       )
        group.addParam('doCC3d', BooleanParam,
                       label='Do 3D radial compensation?',
                       default=True,
                       help='If "No" is selected, the normalized cross correlation (NCC) is made in 2D, '
                            'otherwise radial average is compensated for doing a NCC in 3D.')

        group = form.addGroup('Radial averages',
                              condition='%s or (%s and expertLevel == %s)' %
                                        (AP_CONDITION, NOT_AP_CONDITION, LEVEL_ADVANCED))
        group.addParam('ccMetric', EnumParam,
                       label='Cross correlation metric',
                       choices=['Cross-correlation within the mask',
                                'Mask normalized similarity',
                                'Full cross-correlation'],
                       default=CC_WITHIN_MASK,
                       condition=AP_CONDITION,
                       help='Metric used when computing the cross correlation matrix among the 2D particles. '
                            'Considerations:\n'
                            '\t- Mask normalized similarity is referred to negative squared Euclidean distance.\n'
                            '\t- Full cross-correlation (slower than ross-correlation within the mask but allows '
                            'small disalignments between particles).')
        group.addParam('pcaComps', IntParam,
                       label='PCA components for dim. reduction',
                       validators=[GE(0)],
                       default=0,
                       expertLevel=LEVEL_ADVANCED,
                       condition=NOT_AP_CONDITION,
                       help='Number of components (moments) after the reductions.\nIf 0 or None, '
                            'then they will be automatically estimated considering the size ob the input subtomograms.')

        group = form.addGroup('Parameters of chosen classification algorithm')
        group.addParam('aggNClusters', IntParam,
                       label='Number of clusters to find',
                       default=50,
                       condition=NOT_AP_CONDITION)
        group.addParam('apPref', FloatParam,
                       label='Affinity propagation preference (-inf, inf)',
                       default=-6,
                       condition=AP_CONDITION,
                       help='Preference parameter (-inf, inf).\nThe smaller value the higher number of '
                            'potential classes.\nIf None, the median of the affinity class is considered.')
        group.addParam('apDumping', FloatParam,
                       label='Dumping [0.5, 1)',
                       default=0.5,
                       validators=[GE(0.5), LT(1)],
                       expertLevel=LEVEL_ADVANCED,
                       condition=AP_CONDITION,
                       help='Dumping parameter [0.5, 1), it controls convergence speed.')
        group.addParam('apMaxIter', IntParam,
                       label='Maximum number of iterations',
                       default=2000,
                       condition=AP_CONDITION,
                       expertLevel=LEVEL_ADVANCED)
        group.addParam('apConvIter', IntParam,
                       label='Iterations for fitting the convergence criteria',
                       default=40,
                       condition=AP_CONDITION,
                       expertLevel=LEVEL_ADVANCED)
        group.addParam('apReference', EnumParam,
                       label='Reference 2D image used for classes',
                       choices=['Exemplar', 'Average'],
                       default=AVERAGE,
                       condition=AP_CONDITION,
                       expertLevel=LEVEL_ADVANCED)

        group = form.addGroup('Classification post-processing', expertLevel=LEVEL_ADVANCED)
        group.addParam('apPartSizeFilter', IntParam,
                       label='Minimum number of particles per class',
                       validators=[GE(0)],
                       default=0,
                       expertLevel=LEVEL_ADVANCED,
                       help='Purge classes with less than the specified number of particles. '
                            'If 0, this filter will not be applied.')
        group.addParam('apCCRefFilter', FloatParam,
                       label='Cross-correlation against AP reference filter',
                       validators=[GE(0)],
                       default=0,
                       condition=AP_CONDITION,
                       expertLevel=LEVEL_ADVANCED,
                       help='Purge classes with the cross correlation against the reference '
                            'lower than the specified value. If 0, this filter will not be applied.')
        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.pysegPlaneAlignClassification)
        self._insertFunctionStep(self.createOutputStep)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        subtomoSet = self.inputSubtomos.get()
        subTomoStar = self._getExtraPath(self.inStarName)
        writer = createWriterTomo(isPyseg=True)
        writer.subtomograms2Star(subtomoSet, subTomoStar)
        # Convert the mask format if necessary
        checkMaskFormat(self.inMask.get())

    def pysegPlaneAlignClassification(self):
        # Script called
        Plugin.runPySeg(self, PYTHON, self. _getCommand())

    def createOutputStep(self):
        # Read generated star file and create the output objects:
        # 1) Set of subtomograms
        inSubtomoSet = self.inputSubtomos.get()
        outSubtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
        outSubtomoSet.copyInfo(self.inputSubtomos.get())
        warningMsg, self._dataTable = readPysegSubtomograms(self._getGatheredStarFile(),
                                                            inSubtomoSet,
                                                            outSubtomoSet)
        if warningMsg:
            self._warningMsg.set(warningMsg)
            self._store()

        # 2) Set of classes subtomograms
        classesSet = SetOfClassesSubTomograms.create(self._getPath(), template='setOfClasses%s.sqlite')
        classesSet.setImages(outSubtomoSet)
        self._fillClasses(classesSet)

        self._defineOutputs(**{outputObjects.subtomograms.name: outSubtomoSet,
                               outputObjects.classes.name: classesSet})
        self._defineSourceRelation(inSubtomoSet, outSubtomoSet)
        self._defineSourceRelation(inSubtomoSet, classesSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        return [SEE_METHODS_TAB]

    def _methods(self):
        summary = []
        if self.isFinished():
            sizePostPorcessing = self.apPartSizeFilter.get()
            ccPostProcessing = self.apCCRefFilter.get()
            r3dMsg = 'Radial compensation for 3D'
            r3dMsg = r3dMsg if self.doCC3d.get() else 'No ' + r3dMsg
            summary.append(
                '\n*Particles pre-processing:*\n'
                '   - Low pass Gaussian filter sigma: %i voxels.\n'
                '   - %s' % (self.filterSize.get(), r3dMsg)
            )

            if self.clusteringAlg.get() == AFFINITY_PROP:
                msg = '   - Cross correlation metric: %s\n' % self._decodeCCMetric()
            else:
                msg = '   - Number of components for PCA dimensionality reduction: %i\n' % self.pcaComps.get()
            summary.append(
                '\n*Distance metric calculation:*\n'
                '   - Distance metric: %s\n%s' % (self._decodeDistanceMetric(), msg))

            summary.append(
                '*Classification:*\n'
                '   - Clustering algorithm: %s\n' % self._decodeClusteringAlg()
            )
            if sizePostPorcessing or ccPostProcessing:
                msg = ''
                if sizePostPorcessing:
                    msg += '   - Classes containing less than %i particles were purged\n' % sizePostPorcessing
                if ccPostProcessing:
                    msg += '   - Classes with the cross correlation against the reference lower than %1.2f were ' \
                           'purged' % ccPostProcessing
                summary.append('Post-processing:\n%s' % msg)
        return summary

    def _validate(self):
        errors = []
        tol = 1e-3
        inMask = self.inMask.get()
        inSubtomos = self.inputSubtomos.get()
        maskRes = inMask.getSamplingRate()
        subTomosRes = inSubtomos.getSamplingRate()
        nSubtomos = inSubtomos.getSize()
        xs, ys, zs = inSubtomos.getDimensions()
        xm, ym, zm = inMask.getDimensions()
        if abs(maskRes - subTomosRes) > tol:
            errors.append('Sampling rate of the input subtomograms and the input mask should be the same\n'
                          '%2.3f != %2.3f' % (subTomosRes, maskRes))
        if (xs, ys, zs) != (xm, ym, zm):
            errors.append('The dimensions of the subtomograms and the mask introduced must be the same:\n'
                          '\t- Subtomograms: (x, y, z) = (%i, %i, %i)\n'
                          '\t- Mask: (x, y, z) = (%i, %i, %i)\n' % (xs, ys, zs, xm, ym, zm))
        if self.clusteringAlg.get() != AFFINITY_PROP:
            if self.aggNClusters.get() <= 0 or self.aggNClusters.get() > nSubtomos:
                errors.append('Number of clusters to find must be in range (0, nParticles).')
            if self.pcaComps.get() > nSubtomos:
                errors.append('Number of PCA components must be between 0 and min(n_samples, n_features).')
        return errors

    # --------------------------- UTIL functions -----------------------------------

    def _initialize(self):
        self._outDir = self._getExtraPath(PLANE_ALIGN_CLASS_OUT)
        # Generate output subtomo dir
        makePath(self._outDir)

    def _getGatheredStarFile(self):
        return glob.glob(join(self._outDir, '*_gather.star'))[0]

    def _getNumberOfClasses(self):
        return len(glob.glob(join(self._outDir, '*_split.star')))

    def _getCommand(self):
        alg = self.clusteringAlg.get()
        classCmd = ' '
        classCmd += '%s ' % Plugin.getHome(PLANE_ALIGN_CLASS_SCRIPT)
        classCmd += '--inRootDir scipion '
        classCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        classCmd += '--inMask %s ' % getFinalMaskFileName(self.inMask.get())
        classCmd += '--outDir %s ' % self._outDir
        classCmd += '--filterSize %s ' % self.filterSize.get()
        classCmd += '--procLevel %s ' % (FULL_CLASSIFICATION + 1)  # Numbered from 1 in pyseg
        classCmd += '--doCC3d %s ' % self.doCC3d.get()
        classCmd += '--ccMetric %s ' % self._decodeCCMetric()
        classCmd += '--clusteringAlg %s ' % self._decodeClusteringAlg()
        classCmd += '--distanceMetric %s ' % self._decodeDistanceMetric()
        if alg == AFFINITY_PROP:
            classCmd += '--apPref %s ' % self.apPref.get()
            classCmd += '--apDumping %s ' % self.apDumping.get()
            classCmd += '--apMaxIter %s ' % self.apMaxIter.get()
            classCmd += '--apConvIter %s ' % self.apConvIter.get()
            classCmd += '--apCCRefFilter %s ' % self.apCCRefFilter.get()
        else:
            classCmd += '--pcaComps %s ' % self._estimatePCAComps()
            if alg == AGGLOMERATIVE:
                classCmd += '--aggNClusters %s ' % self.aggNClusters.get()
            else:
                classCmd += '--kmeansNClusters %s ' % self.aggNClusters.get()

        classCmd += '--apPartSizeFilter %s ' % self.apPartSizeFilter.get()
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
        # It deopends on the clustering algorithm, being CC if AP and vectors otherwise
        distanceMetric = self.clusteringAlg.get()
        if distanceMetric == AFFINITY_PROP:
            res = 'ncc_2dz'
        else:
            res = 'vectors'

        return res

    def _fillClasses(self, classesSet):
        classesSet.classifyItems(updateItemCallback=self._updateParticle,
                                 updateClassCallback=self._updateClass,
                                 itemDataIterator=self._dataTable.__iter__())

    @staticmethod
    def _updateParticle(item, row):
        item.setClassId(int(row.rlnClassNumber)+1)  # pyseg classes are 0-indexed

    def _updateClass(self, item):
        classId = item.getObjId()
        fn = self._getReferenceImage(classId)
        item.setAlignment3D()
        item.getRepresentative().setLocation(fn + ':mrc')

    def _getReferenceImage(self, classId):
        # Search in exemplars directory. If file does not exist, search in averages directory
        classId -= 1  # 0-based
        representativesLocation = None
        exemplarssLocation = glob.glob(join(self._outDir, '*_exemplars'))
        averagesLocation = glob.glob(join(self._outDir, '*_averages'))
        if exemplarssLocation:
            if exists(abspath(exemplarssLocation[0])):
                representativesLocation = exemplarssLocation[0]
        else:
            if exists(abspath(averagesLocation[0])):
                representativesLocation = averagesLocation[0]

        return glob.glob(join(representativesLocation, 'class_k%i.mrc' % classId))[0]

    def _estimatePCAComps(self):
        pcaComps = self.pcaComps.get()
        if not pcaComps:
            x, y, _ = self.inputSubtomos.get().getDimensions()
            pcaComps = round(1/100 * 0.5 * x * y)

        return pcaComps

