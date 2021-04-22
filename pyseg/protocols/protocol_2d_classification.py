from pwem.protocols import EMProtocol, PointerParam
from pyworkflow import BETA
from pyworkflow.protocol import String, EnumParam, IntParam, LEVEL_ADVANCED, GT, FloatParam, GE, LT, \
    BooleanParam
from pyworkflow.utils import Message, makePath
from reliontomo.convert import writeSetOfSubtomograms
from scipion.constants import PYTHON
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_SCRIPT, POST_REC_OUT
from pyseg.convert import readStarFile, RELION_SUBTOMO_STAR

# Processing level choices
PARTICLE_FLATENNING = 0  # Particle flattening
CC_MATRIX = 1            # Cross correlation matrix / Feature vectors
FULL_CLASSIFICATION = 2  # Full classification

# Cross correlation metric choices
CC_WITHIN_MASK = 0  # Cross-correlation within mask
SIMILARITY = 1      # Mask normalized similarity (negative squared Euclidean distance)
FULL_CC = 2         # Full cross-correlation

# Clustering algorithm choices
AFFINITY_PROP = 0   # Affinity propagation
AGGLOMERATIVE = 1   # Agglomerative clustering
KMEANS = 2          # K-means

# Affinity propagation condition
AP_CONDITION = 'clusteringAlg == %s and procLevel == %s' % (AFFINITY_PROP, FULL_CLASSIFICATION)

# Affinity propagation distance metric choices
CC_DISTANCE = 0         # Cross correlation used as distance
EUCLIDEAN_DISTANCE = 1  # Euclidean distance among image vectors used as distance metric.

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
        # TODO: Ask Antonio about these parameters to check which are the most critical/important

        group = form.addGroup('Radial averages', condition='procLevel >= %s' % CC_MATRIX)
        group.addParam('ccMetric', EnumParam,
                       label='Cross correlation metric',
                       choices=['Cross-correlation within the mask',
                                'Mask normalized similarity',
                                'Full cross-correlation'],
                       default=CC_WITHIN_MASK,
                       help='Metric used when computing the cross correlation matrix among the 2D particles. '
                            'Considerations:\n'
                            '\t- Mask normalized similarity is referred to negative squared Euclidean distance.\n'
                            '\t- Full cross-correlation (slower than ross-correlation within the mask but allows '
                            'small disalignments between particles).')

        group = form.addGroup('Classification', condition='procLevel >= %s' % FULL_CLASSIFICATION)
        group.addParam('clusteringAlg', EnumParam,
                       label='Clustering algorithm',
                       choices=['Affinity propagation',
                                'Agglomerative clustering',
                                'K-means'],
                       default=AFFINITY_PROP)
        group.addParam('distanceMetric', EnumParam,
                       label='Distance metric',
                       condition='clusteringAlg == %s' % AFFINITY_PROP,
                       choices=['Cross correlation',
                                'Euclidean distance among image vectors'],
                       default=CC_DISTANCE)
        # TODO: in initialization, set distanceMetric to 1 (vectors) if clustering algorithm != AP
        group.addParam('pcaComps', IntParam,
                       label='PCA components for dim. reduction',
                       validators=[GT(0)],
                       default=0,
                       help='Number of components (moments) after the reductions.\nIf None, '
                            'then n_comp == min(n_samples, n_features) - 1')
        group.addParam('aggNClusters', IntParam,
                       label='Number of clusters to find',
                       default=50,
                       condition='clusteringAlg in [%s, %s]' % (AGGLOMERATIVE, KMEANS))
        # TODO: in validation, check that pcaComps > 0 if the clustering algorithm is AGG or K-means

        group = form.addGroup('Affinity propagation',
                              condition=AP_CONDITION)
        group.addParam('apPref', FloatParam,
                       label='Preference (-inf, inf)',
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
                       label='reference 2D image used for classes',
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
                       help='Purging classes with the cross correlation against the reference '
                            'lower than the specified value.')


        # group.addParam('aggLinkage', EnumParam,
        #                label='linkage criterion',
        #                choices=[])
        #TODO: ask Antonio about the linkage criterion, if other values apart from ward are implemented (it seems not to be the case)

        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        outStar = self._getExtraPath(POST_REC_OUT + '.star')
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('pysegPostRec', outStar)
        self._insertFunctionStep('createOutputStep', outStar)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        imgSet = self.inputSubtomos.get()
        imgStar = self._getExtraPath(self.inStarName)
        writeSetOfSubtomograms(imgSet, imgStar, isPysegPosRec=True)

    def pysegPostRec(self, outStar):
        # Generate output subtomo dir
        outDir = self._getExtraPath(POST_REC_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self. _getCommand(outDir, outStar))

    def createOutputStep(self, outStar):
        # Read generated star file and create the output objects
        self.subtomoSet = self._createSetOfSubTomograms()
        self.subtomoSet.setSamplingRate(self.inMask.get().getSamplingRate())
        warningMsg = readStarFile(self, self.subtomoSet, RELION_SUBTOMO_STAR, starFile=outStar)
        if warningMsg:
            self.warningMsg = String(warningMsg)
            self._store()

        self._defineOutputs(outputSetOfSubtomogram=self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('Generated files location:\n'
                           'Subtomograms files directory: %s\n'
                           'Star file: %s' %
                           (self._getExtraPath(POST_REC_OUT),
                            self._getExtraPath(POST_REC_OUT + '.star')))
        return summary

    # --------------------------- UTIL functions -----------------------------------

    def _getCommand(self, outDir, outStar):
        posRecCmd = ' '
        posRecCmd += '%s ' % Plugin.getHome(POST_REC_SCRIPT)
        posRecCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        posRecCmd += '--inMask %s ' % self.inMask.get().getFileName()
        posRecCmd += '--outDir %s ' % outDir
        posRecCmd += '--outStar %s ' % outStar
        posRecCmd += '-j %s ' % self.numberOfThreads.get()
        return posRecCmd
