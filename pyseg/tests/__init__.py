from pyworkflow.tests import DataSet

DataSet(name='emd_10439', folder='emd_10439',
        files={
               'tomoEmd10439': 'tomograms/emd_10439.mrc',
               'coords3dStarFile': 'importFromStarFiles/picking_001_parts.star',
               'coords3dStarFileWithSRate': 'importFromStarFiles/picking_001_parts_with_sRate.star',
               'subtomogramsStarFile': 'importFromStarFiles/class_ap_r_ali_k1_split.star',
               'scipionSqlite3dCoords': 'importFromScipionSqlite/coordinates.sqlite',
               'tomomaskAnnotated': 'tomomasksAnnotated/emd_10439_materials.mrc'
        })

