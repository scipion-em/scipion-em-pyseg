from os import remove
from os.path import join, abspath, exists
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.config import Config
from pyworkflow.utils import magentaStr
from relion.convert import Table
from pyseg.constants import TOMOGRAM, VESICLE, PYSEG_LABEL, MASK
from pyseg.protocols import *


class TestPysegFromPresegToPicking(BaseTest):
    """ """
    ds = DataSet.getDataSet('pyseg')
    preSegPath = join(Config.SCIPION_TESTS, 'pyseg', 'preseg')
    preSegStar = join(preSegPath, 'preseg.star')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        # cls.deltaAng = 1e-4
        # cls.deltaShift = 1e-4

    @classmethod
    def tearDownClass(cls):
        if exists(cls.preSegStar):
            remove(cls.preSegStar)

    def _genPreSegStar(self):
        """Required because pyseg's preseg requires the absolute paths contained
        in the star file"""

        # Headers for Relion's CTS star file
        preSegTable = Table(columns=[TOMOGRAM,
                                     VESICLE,
                                     PYSEG_LABEL,
                                     MASK])

        # Write star file fot 1 tomogram with 3 vesicles
        nVesicles = 3
        tomoName = abspath(self.ds.getFile('presegTomo'))
        mask = abspath(self.ds.getFile('presegMask'))
        for i in range(nVesicles):
            preSegTable.addRow(tomoName,
                               tomoName,
                               i + 1,
                               mask)

        preSegTable.write(self.preSegStar)

    def _runPreseg(self):
        print(magentaStr("\n==> Running preSeg:"))
        self._genPreSegStar()
        protPreseg = self.newProtocol(
            ProtPySegPreSegParticles,
            inStar=self.preSegStar,
            spOffVoxels=30,
            sgVoxelSize=6.87,
            sgMembThk=40,
            sgMembNeigh=200
        )

        protPreseg.setObjLabel('Preseg')
        protPreseg = self.launchProtocol(protPreseg)
        return protPreseg

    def _runGraphs(self, presegProt):
        print(magentaStr("\n==> Running graphs:"))
        protGraphs = self.newProtocol(
            ProtPySegPreSegParticles,
            inSegProt=presegProt,
            pixelSize=6.87,
        )

        protGraphs.setObjLabel('Graphs')
        protGraphs = self.launchProtocol(protGraphs)
        return protGraphs

    def test_workflow(self):
        protPreseg = self._runPreseg()
        protGraphs = self._runGraphs(protPreseg)
