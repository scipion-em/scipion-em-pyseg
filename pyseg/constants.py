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
from os.path import join

PYSEG = 'pySeg'
PYSEG_HOME = 'PYSEG_HOME'
# DEFAULT_VERSION = '1.1.0'
# PYSEG_SOURCE_URL = 'https://github.com/anmartinezs/pyseg_system/archive/refs/tags/v%s.tar.gz' % DEFAULT_VERSION
DEFAULT_VERSION = 'jj_py3'
PYSEG_SOURCE_URL = 'https://github.com/anmartinezs/pyseg_system/archive/refs/heads/%s.zip' % DEFAULT_VERSION

PYSEG_ENV_NAME = 'pySeg-%s' % DEFAULT_VERSION
PYSEG_ENV_ACTIVATION = 'PYSEG_ENV_ACTIVATION'
DEFAULT_ACTIVATION_CMD = 'conda activate %s' % PYSEG_ENV_NAME

# Required files location in pyseg-system
PYSEG_SYSTEM_MAIN = 'pyseg_system-%s' % DEFAULT_VERSION
DATA_TUTORIALS = join(PYSEG_SYSTEM_MAIN, 'data', 'tutorials')
FILS_XML_FILES_PATH = join(DATA_TUTORIALS, 'synth_sumb', 'fils', 'in')
CODE_TUTORIALS = join(PYSEG_SYSTEM_MAIN, 'code', 'tutorials')
EXP_SOMB = join(CODE_TUTORIALS, 'exp_somb')
SYNTH_SUMB = join(CODE_TUTORIALS, 'synth_sumb')

PRESEG_SCRIPT = join(EXP_SOMB, 'pre_tomos_seg.py')
GRAPHS_SCRIPT = join(SYNTH_SUMB, 'tracing', 'mb_graph_mp.py')
FILS_SCRIPT = join(SYNTH_SUMB, 'tracing', 'mb_fils_network.py')
FILS_SOURCES = join(FILS_XML_FILES_PATH, 'mb_sources.xml')
FILS_TARGETS = join(FILS_XML_FILES_PATH, 'no_mb_targets.xml')
PICKING_SCRIPT = join(EXP_SOMB, 'mbo_picking.py')
PICKING_SLICES = join(DATA_TUTORIALS, 'exp_somb', 'mb_ext.xml')
POST_REC_SCRIPT = join(SYNTH_SUMB, 'rln', 'post_rec_particles.py')
PLANE_ALIGN_CLASS_SCRIPT = join(SYNTH_SUMB, 'class', 'plane_align_class.py')

# Generated data
GRAPHS_OUT = 'graphs'
FILS_OUT = 'fils'
PICKING_OUT = 'picking'
POST_REC_OUT = 'subtomos_post_rec'
PLANE_ALIGN_CLASS_OUT = 'plane_align_classification'

# Third parties software
CFITSIO = 'cfitsio'
DISPERSE = 'disperse'

# Dataata source codification
FROM_SCIPION = 0
FROM_STAR_FILE = 1

# Star file fields #####################################################################################################
NOT_FOUND = 'Not found'

# Preseg_pre_centered
TOMOGRAM = 'rlnMicrographName'
VESICLE = 'rlnImageName'
MASK = 'psSegImage'
PYSEG_ROT = 'psSegRot'
PYSEG_TILT = 'psSegTilt'
PYSEG_PSI = 'psSegPsi'
PYSEG_OFFSET_X = 'psSegOffX'
PYSEG_OFFSET_Y = 'psSegOffY'
PYSEG_OFFSET_Z = 'psSegOffZ'

# Preseg_centered
PYSEG_LABEL = 'psSegLabel'
RLN_ORIGIN_X = 'rlnOriginX'
RLN_ORIGIN_Y = 'rlnOriginY'
RLN_ORIGIN_Z = 'rlnOriginZ'
