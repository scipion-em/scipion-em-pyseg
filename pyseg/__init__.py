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
    PYSEG_ENV_NAME

_logo = "icon.png"
_references = ['AMartinez-Sanchez2020']
__version__ = '3.0.0'


class Plugin(pwem.Plugin):

    _homeVar = PYSEG_HOME

    @classmethod
    def _defineVariables(cls):
        # cryoCARE does NOT need EmVar because it uses a conda environment.
        cls._defineVar(PYSEG_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)

    @classmethod
    def getPysegEnvActivation(cls):
        activation = cls.getVar(PYSEG_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep
        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch pyseg. """
        environ = Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']

        return environ

    @classmethod
    def defineBinaries(cls, env):
        PYSEG_INSTALLED = '%s_installed' % PYSEG
        scipion_home = os.environ.get("SCIPION_HOME", None)
        em_root = os.environ.get("EM_ROOT", None)
        zipFile = join(scipion_home, em_root, PYSEG, basename(PYSEG_SOURCE_URL))

        # Try to get CONDA activation command
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

        # Download and extract PySeg Source code
        installationCmd += ' && wget ' + PYSEG_SOURCE_URL
        installationCmd += ' && unzip -q %s' % zipFile
        installationCmd += ' && touch %s' % PYSEG_INSTALLED  # Flag installation finished

        pyseg_commands = [(installationCmd, PYSEG_INSTALLED)]

        env.addPackage(PYSEG,
                       tar='void.tgz',
                       commands=pyseg_commands,
                       neededProgs=["wget"],
                       default=True)

    @classmethod
    def runPySeg(cls, protocol, program, args, cwd=None):
        """ Run pySeg command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(),
                                       cls.getPysegEnvActivation(),
                                       program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)




