import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

class Dedisperser:
    """
    Dedisperses data based on DM and reference frequency.
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
        Adjust MJD based on reference frequency and padding.
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
            raise RuntimeError(f"Geometry not matching for dedispersion: fb.shape = {fb.shape}, expected (T, F) with F = {self.nchans}")

        osamples = nsamples - self.delay_umax
        db = np.zeros((osamples, self.nchans), dtype=np.float32)

        for ichan, delay in enumerate(self.delay_u):
            ib = delay
            jb = ib + osamples
            db[..., ichan] = fb[ib:jb, ichan]

        return db

def plot_dynamic_spectrum(fits_file, dm, ref_freq, output_png=None):
    """
    Plots the stitched summed dynamic spectrum (Stokes 1 + Stokes 2) for a given FITS file,
    with DAT_SCL and DAT_OFFS applied, and performs dedispersion.

    Parameters:
    - fits_file: str, path to the input FITS file.
    - dm: float, dispersion measure to be used for dedispersion.
    - ref_freq: float, reference frequency (MHz) for dedispersion.
    - output_png: str, path to save the output plot as a PNG file (optional). If None, the plot is shown interactively.
    """
    with fits.open(fits_file) as hdulist:
        # Access the binary table containing the 'DAT_SCL' and 'DAT_OFFS' columns
        data_table = hdulist[1].data  # Typically, the binary table is in the second HDU (index 1)

        # Extract the 'DAT_SCL' and 'DAT_OFFS' columns
        dat_scl = data_table['DAT_SCL']
        dat_scl = dat_scl.reshape((5, 1280, 4))  # Reshape to (subint, stokes, frequency channels)
        dat_offs = data_table['DAT_OFFS']
        dat_offs = dat_offs.reshape((5, 1280, 4))  # Reshape to (subint, stokes, frequency channels)

        # Initialize an empty list to hold the corrected dynamic spectra for each subintegration
        all_summed_spectra = []

        # Loop through each subintegration (assuming there are 5 subintegrations)
        for subint in range(5):
            # Extract the dynamic spectrum for the first two Stokes parameters
            dynamic_spectrum_stokes1 = hdulist[1].data['DATA'][subint, :, 3, :]  # First Stokes parameter (time, freq) AA
            dynamic_spectrum_stokes2 = hdulist[1].data['DATA'][subint, :, 3, :]  # Second Stokes parameter (time, freq) BB

            # Apply DAT_SCL and DAT_OFFS to both Stokes parameters
            corrected_spectrum_stokes1 = (dynamic_spectrum_stokes1 - 127.5) * dat_scl[subint, :, 3] + dat_offs[subint, :, 3]
            corrected_spectrum_stokes2 = (dynamic_spectrum_stokes2 - 127.5) * dat_scl[subint, :, 3] + dat_offs[subint, :, 3]

            # Sum the corrected Stokes parameters
            summed_dynamic_spectrum = -corrected_spectrum_stokes1 -corrected_spectrum_stokes2

            # Append the summed dynamic spectrum to the list
            all_summed_spectra.append(summed_dynamic_spectrum)

        # Concatenate all subintegrations along the time axis
        stitched_summed_spectrum = np.concatenate(all_summed_spectra, axis=0)

        # Normalize the dynamic spectrum
        stitched_summed_spectrum = (stitched_summed_spectrum - np.mean(stitched_summed_spectrum, axis=0)) / np.std(stitched_summed_spectrum, axis=0)

        # Frequency information from the FITS file
        freqs = hdulist[1].data['DAT_FREQ']
        print("Frequency shape:", freqs.shape)  # Print the shape to understand its structure

        if freqs.ndim == 2:
            freqs = freqs[0, :]  # Adjust based on your specific FITS file structure
        elif freqs.ndim == 1:
            freqs = freqs  # Already in the expected 1D form
        else:
            raise ValueError("Unexpected shape for frequency data")

    # Dedisperser initialization
    tsamp = hdulist[1].header['TBIN']  # Sampling time in seconds
    dedisperser = Dedisperser(dm, tsamp, freqs, ref_freq=ref_freq)
    
    # Apply dedispersion
    dedispersed_spectrum = dedisperser(stitched_summed_spectrum)

    # Plot the stitched summed dynamic spectrum
    print(dedispersed_spectrum.shape)
    time_bins = np.arange(dedispersed_spectrum.shape[0])  # Total time bins after stitching
    frequency_channels = np.arange(dedispersed_spectrum.shape[1])  # Frequency channels

    plt.figure(figsize=(10, 8))
    plt.imshow(dedispersed_spectrum.T, aspect='auto', origin='lower', cmap='plasma',
               extent=[time_bins.min(), time_bins.max(), frequency_channels.min(), frequency_channels.max()],
               vmin=np.nanpercentile(dedispersed_spectrum.T, 5), vmax=np.nanpercentile(dedispersed_spectrum.T, 95))
    plt.colorbar(label='Intensity')
    plt.xlabel('Time Bins')
    plt.ylabel('Frequency Channels')
    plt.title('Total Intensity Dynamic Spectrum')
    plt.tight_layout()

    if output_png:
        plt.savefig(output_png)
        plt.close()
    else:
        plt.show()

# Example usage
if __name__ == "__main__":
    plot_dynamic_spectrum('b0355_band1.fits', dm=57, ref_freq=2100, output_png='b0355_test.png')

