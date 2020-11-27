# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

PYSEG = 'pySeg'
PYSEG_HOME = 'PYSEG_HOME'
BRANCH = 'nightly'
PYSEG_SOURCE_URL = 'https://github.com/anmartinezs/pyseg_system/archive/%s.zip' % BRANCH

PYSEG_ENV_NAME = 'pySeg_env'
PYSEG_ENV_ACTIVATION = 'PYSEG_ENV_ACTIVATION'
DEFAULT_ACTIVATION_CMD = 'conda activate %s' % PYSEG_ENV_NAME

# Generated data
POST_REC_OUT = 'subtomos_post_rec'

# Third parties software
CFITSIO = 'cfitsio'
DISPERSE = 'disperse'