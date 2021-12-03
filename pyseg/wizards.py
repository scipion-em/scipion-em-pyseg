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
from pwem.wizards import EmWizard
from pyseg.protocols import ProtPySegGraphs


class PysegGraphsDistanceToMbWizard(EmWizard):
    maxLenParamName = 'maxLen'
    _targets = [(ProtPySegGraphs, [maxLenParamName])]

    def show(self, form):
        pysegGraphsProt = form.protocol
        presegProt = pysegGraphsProt.inSegProt.get()
        if not presegProt:
            print('An input pre-segmentation protocol is required to get the automated value for '
                  'the current parameter.')
            return

        form.setVar(self.maxLenParamName, presegProt.sgMembNeigh.get())
