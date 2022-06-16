# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pyworkflow.tests import DataSet

EMD_10439 = 'emd_10439'


class DataSetEmd10439(Enum):
    tomoEmd10439 = 'tomograms/emd_10439.mrc'
    coords3dStarFile = 'importFromStarFiles/picking_001_parts.star'
    coords3dStarFileWithSRate = 'importFromStarFiles/picking_001_parts_with_sRate.star'
    subtomogramsStarFile = 'importFromStarFiles/class_ap_r_ali_k1_split.star'
    scipionSqlite3dCoords = 'importFromScipionSqlite/coordinates.sqlite'
    tomoMaskAnnotated = 'tomomasksAnnotated/emd_10439_materials.mrc'


DataSet(name=EMD_10439, folder=EMD_10439, files={el.name: el.value for el in DataSetEmd10439})
