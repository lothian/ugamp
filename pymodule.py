import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def run_ugamp(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    ugacc can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('ugamp')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_global_option('WFN', 'CCSD')
    scf_helper(name, **kwargs)
    psi4.transqt2()
    if ('wfn' in kwargs):
        if (kwargs['wfn'] == 'ccsd'):
            psi4.set_global_option('WFN', 'CCSD')
        elif (kwargs['wfn'] == 'ccsd(t)'):
            psi4.set_global_option('WFN', 'CCSD_T')
        elif (kwargs['wfn'] == 'mp2'):
            psi4.set_global_option('WFN', 'MP2')
        elif (kwargs['wfn'] == 'mp3'):
            psi4.set_global_option('WFN', 'MP3')
        elif (kwargs['wfn'] == 'mp4'):
            psi4.set_global_option('WFN', 'MP4')
    returnvalue = psi4.plugin('ugamp.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)

def run_ugamp_gradient(name, **kwargs):
    psi4.set_global_option('DERTYPE', 'FIRST')
    run_ugamp(name, **kwargs)

# Integration with driver routines
procedures['energy']['ugamp'] = run_ugamp
procedures['gradient']['ugamp'] = run_ugamp_gradient


def exampleFN():
    # Your Python code goes here
    pass
