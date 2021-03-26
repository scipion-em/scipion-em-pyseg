from pyworkflow.tests import DataSet

DataSet(name='pyseg', folder='pyseg',
        files={'presegTomo': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.mrc',
               'presegMask': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_material.mrc',
               'posRecMask': 'mask88.mrc',
               'subtomosTbl': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.tbl',
               'coordsTbl': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.tbl',
               'posRecCoordsDir': 'posRecCoords',
               'coords': 'posRecCoords/*.mrc'
        })
