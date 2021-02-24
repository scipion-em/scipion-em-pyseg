from pyworkflow.tests import DataSet

DataSet(name='pyseg', folder='pyseg',
        files={'presegDir': 'preseg',
               'presegTomo': 'preseg/Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.mrc',
               'presegMask': 'preseg/Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_material.mrc'
        })
