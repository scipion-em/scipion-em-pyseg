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

from pyworkflow.tests import BaseTest, setupTestProject
from pyseg.protocols import *
from tomo.tests import DataSet


class TestPysegImportSubTomograms(BaseTest):
    """ This class check if the protocol to import sub tomograms works
     properly."""
    @classmethod
    def setUpClass(cls):
         setupTestProject(cls)
         cls.dataset = DataSet.getDataSet('tomo-em')
         cls.star = cls.dataset.getFile('pyseg_after_2dcl.star')

    def _runImportPySegSubTomograms(self):

        print(self.star)
        protImport = self.newProtocol(ProtPySegImportSubtomos,
                                      samplingRate=1.35,
                                      starFile=self.star)
        self.launchProtocol(protImport)
        return protImport

    def test_import_pyseg_subtomograms(self):
        protImport = self._runImportPySegSubTomograms()
        output = getattr(protImport, 'outputSubTomograms', None)
        self.assertTrue(output.getSize() == 2)
        self.assertTrue(output.getSamplingRate() == 1.35)
        self.assertTrue(output.getDim()[0] == 128)
        self.assertTrue(output.getDim()[1] == 128)
        self.assertTrue(output.getDim()[2] == 128)
        self.assertTrue(output.getFirstItem().getSamplingRate() == 1.35)
        self.assertTrue(output.getFirstItem().getDim()[0] == 128)
        self.assertTrue(output.getFirstItem().getDim()[1] == 128)
        self.assertTrue(output.getFirstItem().getDim()[2] == 128)
        self.assertTrue(output.getFirstItem().getClassId() == 30)
        self.assertTrue(output.getFirstItem().getCoordinate3D().getX() == 951)
        self.assertTrue(output.getFirstItem().getCoordinate3D().getY() == 800)
        self.assertTrue(output.getFirstItem().getCoordinate3D().getZ() == 478)
