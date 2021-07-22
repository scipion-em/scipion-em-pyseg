from pyworkflow.tests import DataSet

DataSet(name='pyseg', folder='pyseg',
        files={'presegTomo': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.mrc',
               'presegMask': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_material.mrc',
               'posRecMask': 'cylMaskk40.mrc',
               'mbMask': 'mbMask40.mrc',
               'subtomosTbl': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.tbl',
               'coordsTbl': 'Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED.tbl',
               'posRecCoordsDir': 'posRecCoords',
               'coords': 'posRecCoords/*.mrc',
               'tiltSeries': 'Pertuzumab_1_defocus_25um_tomo_7_ali.mrc',
               'tltFile': 'Pertuzumab_1_defocus_25um_tomo_7_ali.tlt',
               'doseFile': 'Pertuzumab_1_defocus_25um_tomo_7_ali_ExpDose.txt'
        })
