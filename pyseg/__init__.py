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

from os.path import join, basename
import pwem
import os

import pyworkflow
from pyworkflow.utils import Environ
from pyseg.constants import PYSEG_HOME, PYSEG, PYSEG_SOURCE_URL, PYSEG_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD, \
    PYSEG_ENV_NAME, CFITSIO, DISPERSE, DEFAULT_VERSION

_logo = "icon.png"
_references = ['MartinezSanchez2020']
__version__ = '3.0.0'


class Plugin(pwem.Plugin):

    _homeVar = PYSEG_HOME

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(PYSEG_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)
        cls._defineEmVar(PYSEG_HOME, PYSEG + '-' + DEFAULT_VERSION)

    @classmethod
    def getPysegEnvActivation(cls):
        activation = cls.getVar(PYSEG_ENV_ACTIVATION)
        scipionHome = pyworkflow.Config.SCIPION_HOME + os.path.sep
        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch pyseg. """
        environ = Environ(os.environ)
        pySegDir = cls.getHome()

        # Add required disperse path to PATH and pyto path to PYTHONPATH
        environ.update({'PATH': join(pySegDir, '%s_build' % DISPERSE, 'bin'),
                        'PYTHONPATH': join(pySegDir, 'pyseg_system-%s' % DEFAULT_VERSION, 'code')
                        })

        return environ

    @classmethod
    def defineBinaries(cls, env):
        # At this point of the installation execution cls.getHome() is None, so the em path should be provided
        pysegHome = join(pwem.Config.EM_ROOT, PYSEG + '-' + DEFAULT_VERSION)

        PYSEG_INSTALLED = '%s_installed' % PYSEG
        thirdPartyPath = join(pysegHome, 'pyseg_system-%s' % DEFAULT_VERSION, 'sys', 'install',
                              DISPERSE, '0.9.24_pyseg_gcc7', 'sources')

        # PySeg Conda environment
        genPySegCondaEnvCmd = cls._genCmdToDefineSegCondaEnv()
        # PySeg source code
        getPySegCmd = cls._genCmdToGetPySegSrcCode(pysegHome)
        # Third party software - CFitsIO
        installCFitsIOCmd, CFITSIO_BUILD_PATH = cls._genCmdToInstallCFitsIO(thirdPartyPath, pysegHome)
        # Third party software - disperse
        installDisperseCmd = cls._genCmdToInstallDisperse(thirdPartyPath, CFITSIO_BUILD_PATH, pysegHome)

        installationCmd = ' && '.join([genPySegCondaEnvCmd, getPySegCmd, installCFitsIOCmd, installDisperseCmd])
        installationCmd += ' && touch %s' % PYSEG_INSTALLED  # Flag installation finished
        env.addPackage(PYSEG,
                       version=DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, PYSEG_INSTALLED)],
                       neededProgs=["wget", "make", "cmake", "tar", "libgsl-dev"],
                       default=True)

    @classmethod
    def _genCmdToGetPySegSrcCode(cls, pySegDir):
        installationCmd = 'wget ' + PYSEG_SOURCE_URL
        installationCmd += ' && unzip %s' % join(pySegDir, basename(PYSEG_SOURCE_URL))
        installationCmd += ' && rm -rf %s' % join(pySegDir, DEFAULT_VERSION + '.zip')  # rm downloaded file (>600 MB)

        return installationCmd

    @classmethod
    def _genCmdToInstallCFitsIO(cls, thirdPartyPath, pySegDir):
        CFITSIO_BUILD_PATH = join(pySegDir, '%s_build' % CFITSIO)
        installationCmd = 'tar zxf %s -C %s && ' % (join(thirdPartyPath, 'cfitsio_3.380.tar.gz'), pySegDir)
        installationCmd += 'cd %s && ' % join(pySegDir, CFITSIO)
        installationCmd += 'mkdir %s && ' % CFITSIO_BUILD_PATH
        installationCmd += './configure --prefix=%s && ' % CFITSIO_BUILD_PATH
        installationCmd += 'make && make install'

        return installationCmd, CFITSIO_BUILD_PATH

    @classmethod
    def _genCmdToInstallDisperse(cls, thirdPartyPath, CFITSIO_BUILD_PATH, pySegDir):
        DISPERSE_BUILD_PATH = join(pySegDir, '%s_build' % DISPERSE)
        installationCmd = 'tar zxf %s -C %s && ' % \
                          (join(thirdPartyPath, 'disperse_v0.9.24_pyseg_gcc7.tar.gz'), pySegDir)
        installationCmd += 'cd %s && ' % join(pySegDir, DISPERSE)
        installationCmd += 'mkdir %s && ' % DISPERSE_BUILD_PATH
        installationCmd += 'cmake . -DCMAKE_INSTALL_PREFIX=%s -DCFITSIO_DIR=%s && ' \
                           % (DISPERSE_BUILD_PATH, CFITSIO_BUILD_PATH)
        installationCmd += 'make && make install && '
        installationCmd += 'cd ..'

        return installationCmd

    @classmethod
    def _genCmdToDefineSegCondaEnv(cls):
        installationCmd = cls.getCondaActivationCmd()

        # Create the environment
        installationCmd += 'conda create -y -n %s  -c conda-forge -c anaconda python=3.6 ' \
                           'opencv=4.2.0 ' \
                           'graph-tool=2.29 && ' % PYSEG_ENV_NAME

        # Activate new the environment
        installationCmd += 'conda activate %s && ' % PYSEG_ENV_NAME

        # Install non-conda required packages
        installationCmd += 'pip install beautifulsoup4==4.9.3 && '
        installationCmd += 'pip install lxml==4.6.3 && '
        installationCmd += 'pip install pillow==6.2.2 &&'
        installationCmd += 'pip install pyfits==3.5 && '
        installationCmd += 'pip install scikit-image==0.14.5 && '
        installationCmd += 'pip install scikit-learn==0.20.4 && '
        installationCmd += 'pip install scikit-fmm==2021.2.2 && '
        installationCmd += 'pip install scipy==1.2.1 && '
        installationCmd += 'pip install vtk==8.1.2 '
        installationCmd += 'pip install astropy==4.1 '
        installationCmd += 'pip install imageio==2.9.0 '
        installationCmd += 'pip install future'

        return installationCmd

    @classmethod
    def runPySeg(cls, protocol, program, args, cwd=None):
        """ Run pySeg command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(),
                                       cls.getPysegEnvActivation(),
                                       program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)





