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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, basename
import pwem
import os
from pyworkflow.utils import Environ
from pyworkflow import Config
from pyseg.constants import PYSEG_HOME, PYSEG, PYSEG_SOURCE_URL, PYSEG_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD, \
    PYSEG_ENV_NAME, CFITSIO, DISPERSE, BRANCH

_logo = "icon.png"
_references = ['AMartinez-Sanchez2020']
__version__ = '3.0.0'


class Plugin(pwem.Plugin):

    _homeVar = PYSEG_HOME

    @classmethod
    def _defineVariables(cls):
        # PySeg does NOT need EmVar because it uses a conda environment.
        cls._defineVar(PYSEG_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)
        cls._defineVar(PYSEG_HOME, PYSEG)

    @classmethod
    def getPysegEnvActivation(cls):
        activation = cls.getVar(PYSEG_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep
        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch pyseg. """
        environ = Environ(os.environ)
        pySegDir = cls._getPySegDir()

        # Add required disperse path to PATH and pyto path to PYTHONPATH
        environ.update({'PATH': join(pySegDir, '%s_build' % DISPERSE, 'bin'),
                        'PYTHONPATH': join(pySegDir, 'pyseg_system-%s' % BRANCH, 'code')
                        })

        return environ

    @classmethod
    def defineBinaries(cls, env):
        pySegDir = cls._getPySegDir()
        thirdPartyPath = join(pySegDir, 'pyseg_system-%s' % BRANCH, 'sys', 'install',
                              DISPERSE, '0.9.24_pyseg_gcc7', 'sources')

        # PySeg source code
        cls._getPySegSourceCode(env, pySegDir)
        # Third party software - CFitsIO
        CFITSIO_BUILD_PATH = cls._compileCFitsIO(env, pySegDir, thirdPartyPath)
        # Third party software - disperse
        cls._compileDisperse(env, pySegDir, thirdPartyPath, CFITSIO_BUILD_PATH)
        # PySeg Conda environment
        cls._creadtePySegCondaEnv(env)

    @staticmethod
    def _getPySegSourceCode(env, pySegDir):
        PYSEG_DOWNLOADED = join(pySegDir, '%s_downloaded' % PYSEG + '_source')
        installationCmd = 'wget ' + PYSEG_SOURCE_URL
        installationCmd += ' && unzip -q %s' % join(pySegDir, basename(PYSEG_SOURCE_URL))
        installationCmd += ' && rm -rf %s' % join(pySegDir, BRANCH + '.zip')  # Remove the downloaded zip file (>600 MB)
        installationCmd += ' && touch %s' % PYSEG_DOWNLOADED  # Flag download finished
        env.addPackage(PYSEG,
                       tar='void.tgz',
                       commands=[(installationCmd, PYSEG_DOWNLOADED)],
                       neededProgs=["wget", "unzip"],
                       default=True)

    @staticmethod
    def _compileCFitsIO(env, pySegDir, thirdPartyPath):
        CFITSIO_BUILD_PATH = join(pySegDir, '%s_build' % CFITSIO)
        CFITSIO_INSTALLED = join(CFITSIO_BUILD_PATH, '%s_installed' % CFITSIO)
        installationCmd = 'tar zxf %s -C %s && ' % (join(thirdPartyPath, 'cfitsio_3.380.tar.gz'), pySegDir)
        installationCmd += 'cd %s && ' % join(pySegDir, CFITSIO)
        installationCmd += 'mkdir %s && ' % CFITSIO_BUILD_PATH
        installationCmd += './configure --prefix=%s && ' % CFITSIO_BUILD_PATH
        installationCmd += 'make && make install'
        installationCmd += ' && touch %s' % CFITSIO_INSTALLED  # Flag installation finished
        env.addPackage(CFITSIO,
                       tar='void.tgz',
                       commands=[(installationCmd, CFITSIO_INSTALLED)],
                       neededProgs=["make", "tar"],
                       default=True)

        return CFITSIO_BUILD_PATH

    @staticmethod
    def _compileDisperse(env, pySegDir, thirdPartyPath, CFITSIO_BUILD_PATH):
        DISPERSE_BUILD_PATH = join(pySegDir, '%s_build' % DISPERSE)
        DISPERSE_INSTALLED = join(DISPERSE_BUILD_PATH, '%s_installed' % DISPERSE)
        installationCmd = 'tar zxf %s -C %s && ' % \
                          (join(thirdPartyPath, 'disperse_v0.9.24_pyseg_gcc7.tar.gz'), pySegDir)
        installationCmd += 'cd %s && ' % join(pySegDir, DISPERSE)
        installationCmd += 'mkdir %s && ' % DISPERSE_BUILD_PATH
        installationCmd += 'cmake . -DCMAKE_INSTALL_PREFIX=%s -DCFITSIO_DIR=%s && ' \
                           % (DISPERSE_BUILD_PATH, CFITSIO_BUILD_PATH)
        installationCmd += 'make && make install'
        installationCmd += ' && touch %s' % DISPERSE_INSTALLED  # Flag installation finished
        env.addPackage(DISPERSE,
                       tar='void.tgz',
                       commands=[(installationCmd, DISPERSE_INSTALLED)],
                       neededProgs=["make", "tar"],
                       default=True)

    @classmethod
    def _creadtePySegCondaEnv(cls, env):
        PYSEG_INSTALLED = '%s_installed' % PYSEG
        installationCmd = cls.getCondaActivationCmd()

        # Create the environment
        installationCmd += 'conda create -y -n %s -c conda-forge -c anaconda python=2.7 ' \
                           'opencv graph-tool && ' % PYSEG_ENV_NAME

        # Activate new the environment
        installationCmd += 'conda activate %s && ' % PYSEG_ENV_NAME

        # Install non-conda required packages
        installationCmd += 'pip install beautifulsoup4 && '
        installationCmd += 'pip install lxml && '
        installationCmd += 'pip install pillow &&'
        installationCmd += 'pip install pyfits && '
        installationCmd += 'pip install scikit-image && '
        installationCmd += 'pip install scikit-learn && '
        installationCmd += 'pip install scikit-fmm && '
        installationCmd += 'pip install scipy && '
        installationCmd += 'pip install vtk'
        installationCmd += ' && touch %s' % PYSEG_INSTALLED  # Flag installation finished

        pyseg_commands = [(installationCmd, PYSEG_INSTALLED)]

        env.addPackage(PYSEG + '_execution_env',
                       tar='void.tgz',
                       commands=pyseg_commands,
                       neededProgs=["wget", "unzip"],
                       default=True)

    @classmethod
    def runPySeg(cls, protocol, program, args, cwd=None):
        """ Run pySeg command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(),
                                       cls.getPysegEnvActivation(),
                                       program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @staticmethod
    def _getPySegDir():
        scipion_home = os.environ.get("SCIPION_HOME", None)
        em_root = os.environ.get("EM_ROOT", None)
        return join(scipion_home, em_root, PYSEG)





