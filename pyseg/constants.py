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
from os.path import join

PYSEG = 'pySeg'
PYSEG_HOME = 'PYSEG_HOME'
BRANCH = 'nightly'
DEFAULT_VERSION = '1.0.%s' % BRANCH
PYSEG_SOURCE_URL = 'https://github.com/anmartinezs/pyseg_system/archive/%s.zip' % BRANCH

PYSEG_ENV_NAME = 'pySeg-%s' % DEFAULT_VERSION
PYSEG_ENV_ACTIVATION = 'PYSEG_ENV_ACTIVATION'
DEFAULT_ACTIVATION_CMD = 'conda activate %s' % PYSEG_ENV_NAME

# Required files location in pyseg-system
PYSEG_SYSTEM_MAIN = 'pyseg_system-%s' % BRANCH
SYNTH_SUMB = join(PYSEG_SYSTEM_MAIN, 'data', 'tutorials', 'synth_sumb')
GRAPHS_SCRIPT = join(SYNTH_SUMB, 'tracing', 'mb_graph_mp.py')
FILS_SCRIPT = join(SYNTH_SUMB, 'tracing', 'mb_fils_network.py')
FILS_SOURCES = join(SYNTH_SUMB, 'fils', 'in', 'mb_sources.xml')
FILS_TARGETS = join(SYNTH_SUMB, 'fils', 'in', 'no_mb_targets.xml')
PICKING_SCRIPT = join(PYSEG_SYSTEM_MAIN, 'code', 'tutorials', 'exp_somb', 'mbo_picking.py')
PICKING_SLICES = join(SYNTH_SUMB, 'data', 'tutorials', 'exp_somb', 'mb_ext.xml')
POST_REC_SCRIPT = join(SYNTH_SUMB, 'rln', 'post_rec_particles.py')

# Generated data
GRAPHS_OUT = 'graphs'
FILS_OUT = 'fils'
PICKING_OUT = 'picking'
POST_REC_OUT = 'subtomos_post_rec'

# Third parties software
CFITSIO = 'cfitsio'
DISPERSE = 'disperse'
