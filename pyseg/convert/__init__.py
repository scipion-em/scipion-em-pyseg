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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
from emtable import Table
from pyseg.convert.convert import PysegStarReader


def createPysegReader(starFile, **kwargs):
    dataTable = Table()
    dataTable.read(starFile, tableName=None)
    return PysegStarReader(starFile, dataTable, **kwargs)


def readPysegCoordinates(starFile, coordSet, precedentsSet):
    reader = createPysegReader(starFile)
    return reader.starFile2Coords3D(coordSet, precedentsSet)


def readPysegSubtomograms(starFile, inSubtomos, outSubtomos):
    reader = createPysegReader(starFile)
    return reader.starFile2Subtomograms(inSubtomos, outSubtomos)

