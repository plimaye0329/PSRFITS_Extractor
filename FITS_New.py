import numpy as np
import matplotlib.pyplot as plt
import fitsio as fio
from skimage.measure import block_reduce
import os

class Dedisperser:
    """
    Dedisperses data based on DM and reference frequency
    """
    DM_FAC = 4.148741601E3  # Constant for dedispersion

    def __init__(self, dm, tsamp, freqs, ref_freq=0):
        self.tsamp = tsamp
        self.freqs = np.array(freqs).copy()
        self.nchans = self.freqs.shape[0]
        self.freqs_max = self.freqs.max()

        # Reference frequency
        self.ref_freq = ref_freq
        self._flag_ref = False
        if self.ref_freq <= 0:
            self.ref_freq = self.freqs_max
            self._flag_ref = True

        self.dm = dm

        # Calculate delays in seconds
        self.delay_s = self.dm * self.DM_FAC * (self.freqs**-2 - self.ref_freq**-2)

        # Convert delays to time units (in terms of sampling time)
        self.delay_u = np.int32(self.delay_s / self.tsamp)
        if np.any(self.delay_u < 0):
            raise RuntimeError("Negative delay encountered, dedispersion not implemented properly")

        self.delay_umax = np.max(np.abs(self.delay_u))

    def ref_time(self, mjd, padding=0.0):
        """
        Adjust MJD based on reference frequency and padding
        """
        if self._flag_ref:
            return mjd - (padding / 86_400)  # Quick return if ref_freq is already the max

        _corr = self.dm * self.DM_FAC * (self.freqs_max**-2 - self.ref_freq**-2)
        return mjd - ((_corr + padding) / 86_400)

    def __call__(self, fb):
        """
        De-disperse the frequency-time data (fb).
        Input: fb with shape (T, F)
        """
        nsamples, nchans = fb.shape
        if nsamples < self.delay_umax or nchans != self.nchans:
            raise RuntimeError("Geometry not matching for dedispersion")

        osamples = nsamples - self.delay_umax
        db = np.zeros((osamples, self.nchans), dtype=np.float32)

        for ichan, delay in enumerate(self.delay_u):
            ib = delay
            jb = ib + osamples
            db[..., ichan] = fb[ib:jb, ichan]

        return db


def normalize(fb):
    """
    Normalize data (TF layout)
    """
    _mean = np.median(fb, axis=0)
    _std = fb.std(axis=0)
    return np.divide(fb - _mean, _std, where=_std != 0)


def plot_dynamic_spectrum(dynamic_spectrum, freqs, times, dedispersed_spectrum=None):
    """
    Plot the dynamic spectrum and dedispersed spectrum (if available).
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    ax1.imshow(block_reduce(dynamic_spectrum, (2, 4), func=np.mean, cval=0.).T, aspect='auto', cmap='plasma', origin='lower', interpolation='none')
    ax1.set_title("Dynamic Spectrum")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Frequency")

    if dedispersed_spectrum is not None:
        ax2.imshow(block_reduce(dedispersed_spectrum, (2, 4), func=np.mean, cval=0.).T, aspect='auto', cmap='plasma', origin='lower', interpolation='none')
        ax2.set_title("De-dispersed Spectrum")
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Frequency")

    plt.tight_layout()
    plt.savefig('b0355_dedispersed.png')
    plt.show()


def main(fits_file, dm, ref_freq, output_dir):
    # Load the fits file
    primary_header = fio.read_header(fits_file, ext=0)
    subint_header = fio.read_header(fits_file, ext='SUBINT')

    tsamp = subint_header['TBIN']
    nchans = subint_header['NCHAN']
    total_samps = subint_header['NSTOT']

    freqs, = fio.read(fits_file, ext='SUBINT', columns=['DAT_FREQ'], rows=1)[0]
    dynamic_spectrum = np.zeros((total_samps, nchans), dtype=np.float32)

    # Read data from FITS file
    for i in range(subint_header['NAXIS2']):
        data = fio.read(fits_file, ext='SUBINT', columns=['DATA'], rows=i)['DATA']
        dynamic_spectrum[i * subint_header['NSBLK']:(i + 1) * subint_header['NSBLK'], :] = data[0, :, 0, :]

    # Normalize the dynamic spectrum
    normalized_spectrum = normalize(dynamic_spectrum)

    # Dedisperse using Dedisperser class
    dedisperser = Dedisperser(dm, tsamp, freqs, ref_freq=ref_freq)
    dedispersed_spectrum = dedisperser(normalized_spectrum)

    # Create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Plot and save the dynamic spectrum and dedispersed spectrum
    plot_dynamic_spectrum(normalized_spectrum, freqs, np.arange(total_samps) * tsamp, dedispersed_spectrum)

    print("Plot saved successfully in", output_dir)


if __name__ == "__main__":
    fits_file = "b0355_band1.fits"
    dm = 100.0
    ref_freq = 2100.0  # MHz
    output_dir = "output_spectra"

    main(fits_file, dm, ref_freq, output_dir)

