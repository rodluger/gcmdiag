from __future__ import division, print_function, absolute_import, unicode_literals
from . import constants, diag, netcdfcombine, plot, rect, utils
from .diag import *
from .plot import *
from .utils import *
from .netcdfcombine import *
from .rect import *

import os
GCMDIAG_SRC = os.path.dirname(os.path.abspath(__file__))
GCMDIAG_OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'output')