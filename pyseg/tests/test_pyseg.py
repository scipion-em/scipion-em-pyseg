# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
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
from os import remove
from os.path import exists

from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.config import Config, join
from pyworkflow.utils import getParentFolder
from relion.convert import Table

from pyseg.convert import RELION_TOMO_LABELS
from pyseg.protocols import *


class TestPysegImportSubTomograms(BaseTest):
    """ This class check if the protocol to import sub tomograms works
     properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.star = join(Config.SCIPION_TESTS, 'tomo-em', 'starTest.star')
        cls._genDataDicts()

    @classmethod
    def tearDownClass(cls):
        if exists(cls.star):
            remove(cls.star)

    def _runImportPySegSubTomograms(self):

        print(self.star)
        protImport = self.newProtocol(ProtPySegImportSubtomos,
                                      samplingRate=1.35,
                                      starFile=self.star)
        self.launchProtocol(protImport)
        return protImport

    @staticmethod
    def _getKeysStar23():
        return [
            'rlnMicrographName',
            'rlnCtfImage',
            'rlnImageName',
            'rlnCoordinateX',
            'rlnCoordinateY',
            'rlnCoordinateZ',
            'rlnAngleTilt',
            'rlnAngleTiltPrior',
            'rlnAnglePsi',
            'rlnAnglePsiPrior',
            'rlnAngleRot',
            'rlnGroupNumber',
            'rlnOriginX',
            'rlnOriginY',
            'rlnOriginZ',
            'rlnClassNumber',
            'rlnNormCorrection',
            'rlnLogLikeliContribution',
            'rlnMaxValueProbDistribution',
            'rlnNrOfSignificantSamples',
            'rlnRandomSubset',
            'rlnMagnification',
            'rlnDetectorPixelSize'
        ]

    @classmethod
    def _loadTestData(cls):
        path = getParentFolder(cls.star)
        return [['tomo1.mrc', 'tomo2.mrc'],
                  ['wedge1.mrc', 'wedge2.mrc'],
                  [join(path, 'import_particle_000003.mrc'), join(path, 'import_particle_000016.mrc')],
                  [951.053406, 1088.109376],
                  [800.000000, 1082.000000],
                  [478.000000, 432.000000],
                  [86.870899, 48.663964],
                  [87.266841, 58.868786],
                  [136.079624, -57.13563],
                  [126.524687, -47.67025],
                  [17.818407, 5.136058],
                  [13, 5],
                  [9.097985, 10.975485],
                  [0.097985, -3.24202],
                  [1.390485, -0.07202],
                  [30, 30],
                  [1.000000, 1.000000],
                  [3383545, 3387058],
                  [0.383647, 0.506495],
                  [3, 3],
                  [1, 2],
                  [10000, 10000],
                  [4.4, 4.4]
                ]

    @classmethod
    def _genDataDicts(cls):
        keys = cls._getKeysStar23()
        values = cls._loadTestData()
        values1 = [val[0] for val in values]
        values2 = [val[1] for val in values]
        cls.dictList = [{key: val for key, val in zip(keys, values1)},
                        {key: val for key, val in zip(keys, values2)}]


    @classmethod
    def _writeTestStarFile(cls, columns):
        tomoTable = Table(columns=columns)
        for d in cls.dictList:
            vals = [d[key] for key in columns]
            tomoTable.addRow(* vals)

        tomoTable.write(cls.star)

    def test_import_pyseg_subtomograms_23_columns(self):
        self._writeTestStarFile(self._getKeysStar23())  # Write the corresponding star file
        protImport = self._runImportPySegSubTomograms()
        output = getattr(protImport, 'outputSubTomograms', None)
        # Check set attribute
        self._checkSet(output)
        # Check subtomo attributes
        classId = [30, 30]

    def test_import_pyseg_subtomograms_14_columns(self):
        usedKeys = self._writeTestStarFile(RELION_TOMO_LABELS)  # Write the corresponding star file
        protImport = self._runImportPySegSubTomograms()
        subtomoSet = getattr(protImport, 'outputSubTomograms', None)
        # Check set attribute
        self._checkSet(subtomoSet)
        # Check subtomo attributes
        path = getParentFolder(self.star)
        x = [951, 1088]
        y = [800, 1082]
        z = [478, 432]
        filenames = [join(path, 'import_particle_000003.mrc'), join(path, 'import_particle_000016.mrc')]
        volNames = ['tomo1.mrc', 'tomo2.mrc']
        wedges = ['wedge1.mrc', 'wedge2.mrc']
        for i, subtomo in enumerate(subtomoSet):
            self.assertEqual(subtomo.getSamplingRate(), 1.35)
            self.assertEqual(subtomo.getDim()[0], 128)
            self.assertEqual(subtomo.getDim()[1], 128)
            self.assertEqual(subtomo.getDim()[2], 128)
            # self.assertEqual(subtomo.getClassId() == 30)
            self.assertEqual(subtomo.getCoordinate3D().getX(), x[i])
            self.assertEqual(subtomo.getCoordinate3D().getY(), y[i])
            self.assertEqual(subtomo.getCoordinate3D().getZ(), z[i])
            self.assertEqual(subtomo.getFileName(), filenames[i])
            self.assertEqual(subtomo.getVolName(), volNames[i])
            self.assertEqual(subtomo.getCoordinate3D()._3dcftMrcFile, wedges[i])

            # print(subtomo.getCoordinate3D()._3dcftMrcFile)

    def _checkSet(self, output):
        self.assertEqual(output.getSize(), 2)
        self.assertEqual(output.getSamplingRate(), 1.35)
        self.assertEqual(output.getDim()[0], 128)
        self.assertEqual(output.getDim()[1], 128)
        self.assertEqual(output.getDim()[2], 128)
        # self.assertEqual(output.getSamplingRate() == 1.35)
        # self.assertEqual(output.getDim()[0] == 128)
        # self.assertEqual(output.getDim()[1] == 128)
        # self.assertEqual(output.getDim()[2] == 128)
        # self.assertEqual(output.getClassId() == 30)
        # self.assertEqual(output.getCoordinate3D().getX() == 951)
        # self.assertEqual(output.getCoordinate3D().getY() == 800)
        # self.assertEqual(output.getCoordinate3D().getZ() == 478)

        # TODO --> comprobación de los valores de los 2 subtomogramas con los dictioanrios de datos