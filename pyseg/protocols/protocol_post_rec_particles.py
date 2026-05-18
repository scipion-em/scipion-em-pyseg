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

from pwem.protocols import EMProtocol, PointerParam
from pyseg.convert import readPysegSubtomograms
from pyseg.utils import getFinalMaskFileName, checkMaskFormat
from pyworkflow.protocol import String, FloatParam, LE, GE
from pyworkflow.utils import Message, makePath
from reliontomo.convert import createWriterTomo
from scipion.constants import PYTHON
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtTomoBase

from pyseg import Plugin
from pyseg.constants import POST_REC_OUT, POST_REC_SCRIPT_MEMB_ATT, SEE_METHODS_TAB


class outputObjects(Enum):
    subtomograms = SetOfSubTomograms


class ProtPySegPostRecParticles(EMProtocol, ProtTomoBase):
    """
    Post-processes reconstructed subtomograms through membrane suppression
    and particle conditioning operations in preparation for downstream
    subtomogram averaging, classification, or structural interpretation.

    AI Generated:

    Post Reconstruction Particles (ProtPySegPostRecParticles) — User Manual
        Overview

        The Post Reconstruction Particles protocol is designed to refine
        already reconstructed subtomograms before advanced structural
        analysis. In cryo-electron tomography workflows, reconstructed
        particles frequently contain surrounding membrane signal,
        heterogeneous background densities, and orientation-dependent
        artifacts that can interfere with alignment accuracy and reduce
        the interpretability of averages. This protocol provides a
        biologically oriented preprocessing stage intended to suppress
        unwanted membrane contributions while preserving the molecular
        information of interest.

        The protocol is particularly useful in membrane-associated
        systems, including vesicles, organelles, viral envelopes, and
        membrane protein assemblies. In these contexts, the membrane may
        dominate the density distribution and bias subsequent alignment
        or classification procedures. By attenuating membrane signal,
        the protocol improves the visibility of embedded or associated
        macromolecular complexes and facilitates more reliable downstream
        processing.

        Inputs and Biological Context

        The protocol requires a set of reconstructed subtomograms and a
        primary mask defining the region of interest to preserve during
        processing. The subtomograms are assumed to represent particles
        already extracted and reconstructed from tomographic data. The
        mask should encompass the biologically relevant density while
        minimizing unrelated solvent or neighboring structures.

        An optional membrane suppression mask can also be provided. This
        mask identifies membrane regions whose signal should be reduced
        or attenuated. Such suppression is especially valuable when the
        membrane itself is not the target of analysis and would otherwise
        dominate rotational alignment or classification procedures.

        For biological applications, careful mask design is essential.
        The primary mask should include the structural core of the
        particle while excluding irrelevant densities. The membrane
        suppression mask should focus only on membrane regions whose
        contribution is considered undesirable for the intended analysis.
        Excessively aggressive masking may remove biologically meaningful
        information, whereas insufficient masking may leave strong
        membrane artifacts unresolved.

        Membrane Suppression Strategy

        The membrane suppression stage attenuates density contributions
        associated with membrane regions. Instead of completely removing
        membrane information in all cases, the protocol allows gradual
        attenuation, enabling users to balance structural preservation
        with background reduction.

        From a biological perspective, this flexibility is important
        because membranes can contain meaningful contextual information.
        In some studies, complete suppression may be desirable when the
        goal is to focus exclusively on soluble or protruding domains.
        In other situations, partial attenuation preserves orientation
        cues or local structural relationships while still reducing
        alignment bias.

        Lower suppression factors produce stronger attenuation and are
        typically preferred when membrane density overwhelms the particle
        signal. Higher factors preserve more membrane information and may
        be advantageous when membrane geometry contributes to biological
        interpretation.

        Data Consistency and Sampling Considerations

        Accurate processing requires consistency between subtomograms and
        masks. All input volumes should share the same sampling rate and
        spatial scaling to ensure biologically meaningful results.
        Mismatched voxel sizes can introduce distortions, incorrect mask
        placement, or inconsistent attenuation effects.

        In practical cryo-ET workflows, masks are often generated from
        segmentation procedures or external processing pipelines. Before
        running this protocol, users should verify that masks align
        correctly with the subtomogram coordinate system and preserve the
        intended structural regions.

        Outputs and Interpretation

        The protocol produces a new set of processed subtomograms ready
        for downstream structural analysis. The resulting particles
        preserve the original biological identity of the dataset while
        incorporating the requested suppression and conditioning
        operations.

        Biologically, the processed particles are generally better suited
        for subtomogram averaging, classification, and alignment because
        distracting membrane signal is reduced. This often improves the
        interpretability of structural heterogeneity and enhances the
        recovery of molecular features that may otherwise remain obscured.

        The protocol also preserves processing metadata describing the
        operations applied during post reconstruction conditioning. This
        information is useful for reproducibility and for comparing
        alternative preprocessing strategies during iterative refinement
        workflows.

        Practical Recommendations

        In most membrane-associated cryo-ET studies, it is advisable to
        begin with moderate membrane suppression rather than complete
        elimination. Excessive suppression may remove biologically
        relevant densities or introduce artificial discontinuities near
        membrane-contact regions.

        Users should visually inspect representative processed particles
        after execution to confirm that the molecular region of interest
        remains intact and that membrane attenuation behaves as expected.
        Iterative optimization of masks and suppression factors is often
        beneficial in challenging datasets with crowded environments or
        highly curved membranes.

        For datasets involving membrane proteins or assemblies tightly
        associated with lipid bilayers, preserving partial membrane
        information may improve biological interpretation and maintain
        meaningful structural context.

        Final Perspective

        In subtomogram analysis workflows, post reconstruction
        conditioning is an important preparatory step that can strongly
        influence downstream structural interpretation. Careful control
        of membrane suppression and mask definition allows researchers to
        enhance particle quality while preserving the biologically
        relevant information needed for accurate alignment, averaging,
        and classification.
    """

    _label = 'posrec'
    inStarName = 'input_particles.star'
    warningMsg = None
    subtomoSet = None

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
                      help='Select the input subtomograms from the project.')
        form.addParam('inMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask',
                      important=True,
                      allowsNull=False,
                      help='Mask used for the post processing')
        group = form.addGroup('Membrane suppression')
        group.addParam('mbMask', PointerParam,
                       pointerClass='VolumeMask',
                       allowsNull=True,
                       label='Mask for membrane supression')
        group.addParam('mbSupFactor', FloatParam,
                       label='Membrane suppression factor [0, 1]',
                       default=0.5,
                       allowsNull=False,
                       validators=[GE(0), LE(1)],
                       condition='mbMask',
                       help='Value 0 suppress the area corresponding to the suppression mask, while higher values '
                            'up to 1 attenuate it.')
        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        outStar = self._getExtraPath(POST_REC_OUT + '.star')
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.pysegPostRec, outStar)
        self._insertFunctionStep(self.createOutputStep, outStar)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        # Check masks format and convert if necessary
        checkMaskFormat(self.inMask.get())  # Subtomogram mask
        if self.mbMask.get():
            checkMaskFormat(self.mbMask.get())  # Membrane mask for attenuation (optional)
        # Write star from set of subtomograms
        subtomoSet = self.inputSubtomos.get()
        subTomoStar = self._getExtraPath(self.inStarName)
        writer = createWriterTomo(isPyseg=True)
        writer.subtomograms2Star(subtomoSet, subTomoStar)
        # Convert the mask format if necessary
        checkMaskFormat(self.inMask.get())

    def pysegPostRec(self, outStar):
        # Generate output subtomo dir
        outDir = self._getExtraPath(POST_REC_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getCommand(outDir, outStar))

    def createOutputStep(self, outStar):
        self.subtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
        self.subtomoSet.copyInfo(self.inputSubtomos.get())
        # Read generated star file and create the output objects
        warningMsg, _ = readPysegSubtomograms(outStar,
                                              self.inputSubtomos.get(),
                                              self.subtomoSet)
        if warningMsg:
            self.warningMsg = String(warningMsg)
            self._store()

        self._defineOutputs(**{outputObjects.subtomograms.name: self.subtomoSet})
        self._defineSourceRelation(self.inputSubtomos.get(), self.subtomoSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('*Generated files location*:\n'
                           '\t- Subtomograms files directory: %s\n'
                           '\t- Star file: %s\n'
                           '%s' %
                           (self._getExtraPath(POST_REC_OUT),
                            self._getExtraPath(POST_REC_OUT + '.star'),
                            SEE_METHODS_TAB))
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            mbMask = self.mbMask.get()
            if mbMask:
                methods.append('*Membrane suppression applied with*\n'
                               '\t- MembraneMask = %s\n'
                               '\t- SuppressionFactor = %1.2f\n' %
                               (mbMask.getFileName(), self.mbSupFactor.get()))
            if self.doGaussLowPassFilter.get():
                methods.append('*Gaussian low pass filter applied with*\n'
                               '\t- CutOffResolution[nm] = %2.1f\n'
                               '\t- AmplitudeCutOff = %1.2f\n'
                               '\t- AppliedToCTF = %s\n' %
                               (self.cutOffRes.get(), self.ampCutOff.get(), self.filterCTF.get()))

        return methods

    def _validate(self):
        validationMsg = []
        tol = 0.01
        subTomosRes = self.inputSubtomos.get().getSamplingRate()
        inMask = self.inMask.get()
        maskRes = inMask.getSamplingRate()
        mbMask = self.mbMask.get()
        if abs(maskRes - subTomosRes) > tol:
            validationMsg.append('Sampling rate of the input subtomograms and the input mask should be the same\n'
                                 '%2.3f != %2.3f' % (subTomosRes, maskRes))
        if mbMask:
            mbMaskRes = mbMask.getSamplingRate()
            if abs(mbMaskRes - subTomosRes) > tol:
                validationMsg.append('Sampling rate of the input subtomograms and the input membrane suppression mask '
                                     'should be the same\n'
                                     '%2.3f != %2.3f' % (subTomosRes, mbMaskRes))
        return validationMsg

    # --------------------------- UTIL functions -----------------------------------

    def _getCommand(self, outDir, outStar):
        posRecCmd = ' '
        posRecCmd += '%s ' % Plugin.getHome(POST_REC_SCRIPT_MEMB_ATT)
        posRecCmd += '--inStar %s ' % self._getExtraPath(self.inStarName)
        posRecCmd += '--inMask %s ' % getFinalMaskFileName(self.inMask.get())
        posRecCmd += '--inMaskMbSup %s ' % (getFinalMaskFileName(self.mbMask.get()) if self.mbMask.get() else 'None')
        posRecCmd += '--mbSupFactor %s ' % (self.mbSupFactor.get() if self.mbSupFactor.get() else '0')
        posRecCmd += '--doGaussLowPass %s ' % False
        posRecCmd += '--outDir %s ' % outDir
        posRecCmd += '--outStar %s ' % outStar
        posRecCmd += '-j %s ' % self.numberOfThreads.get()
        return posRecCmd

