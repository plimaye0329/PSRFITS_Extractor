import numpy as np
import matplotlib.pyplot as plt
import glob
from matplotlib.pyplot import figure
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
import sigpyproc.readers as readers
import sigpyproc.timeseries as timeseries
import sigpyproc.block as block
from rich.pretty import Pretty
from scipy.signal import find_peaks
from scipy.interpolate import make_interp_spline
import os
import warnings
import math

def sigread(filepath, dm):
    """
    Reads a .fil file and returns the dedispersed data.

    This function reads a .fil file using the FilReader class from the `readers` module.
    It then dedisperses the data using the provided dispersion measure (dm). It then normalizes
    the data by subtracting the median and dividing by the standard deviation.

    Parameters:
    filepath (str): The path to the .fil file to be read.
    dm (float): The dispersion measure to be used for dedispersion.
    block_size (int): The size of the block to be read from the .fil file.

    Returns:
    numpy.ndarray: The dedispersed data from the .fil file.
    """
    fil = readers.PFITSReader(filepath)
    header = fil.header
    block_size = fil.len()
    data = fil.read_block(1, block_size)
    data_dedisp = data.dedisperse(dm)
    arr = data_dedisp.copy()
    #with warnings.catch_warnings():
        #warnings.simplefilter("ignore", category=RuntimeWarning)
        #arr -= np.nanmedian(arr, axis=-1)[..., None]
        #arr /= np.nanstd(arr, axis=-1)[..., None]

    return arr



