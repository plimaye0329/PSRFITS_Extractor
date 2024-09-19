import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

def plot_dedispersed_dynamic_spectrum(fits_file, dm=57.14, output_png=None):
    """
    Plots the stitched summed dynamic spectrum (Stokes 1 + Stokes 2) with dedispersion applied.
    
    Parameters:
    - fits_file: str, path to the input FITS file.
    - dm: float, dispersion measure (DM) in pc cm^-3. Default is 57.14.
    - output_png: str, path to save the output plot as a PNG file (optional). If None, the plot is shown interactively.
    """
    # Constants
    sampling_time = 128e-6  # 128 microseconds in seconds
    k_dm = 4.15e3  # Dispersion constant in ms

    with fits.open(fits_file) as hdulist:
        # Access the binary table containing the 'DAT_SCL' and 'DAT_OFFS' columns
        data_table = hdulist[1].data  # Typically, the binary table is in the second HDU (index 1)

        # Extract the 'DAT_SCL' and 'DAT_OFFS' columns
        dat_scl = data_table['DAT_SCL'].reshape((5, 1280, 4))  # (subint, frequency channels, stokes)
        dat_offs = data_table['DAT_OFFS'].reshape((5, 1280, 4))

        # Extract frequency information from the header (assuming 'OBSFREQ' and 'OBSBW' are available)
        freq = 1600 #hdulist[1].header['OBSFREQ']  # Center frequency in MHz
        bw = 650 # hdulist[1].header['OBSBW']  # Bandwidth in MHz
        nchans = 1280  # Number of frequency channels
        freqs = np.linspace(freq - bw / 2, freq + bw / 2, nchans)

        # Initialize an empty list to hold the corrected dynamic spectra for each subintegration
        all_summed_spectra = []

        # Loop through each subintegration (assuming there are 5 subintegrations)
        for subint in range(5):
            # Extract the dynamic spectrum for the first two Stokes parameters (Stokes I + Stokes Q)
            dynamic_spectrum_stokes1 = hdulist[1].data['DATA'][subint, :, 3, :]  # First Stokes parameter (time, freq)
            dynamic_spectrum_stokes2 = hdulist[1].data['DATA'][subint, :, 3, :]  # Second Stokes parameter (time, freq)

            # Apply DAT_SCL and DAT_OFFS to both Stokes parameters
            corrected_spectrum_stokes1 = (dynamic_spectrum_stokes1 - 127.5) * dat_scl[subint, :, 3] + dat_offs[subint, :, 3]
            corrected_spectrum_stokes2 = (dynamic_spectrum_stokes2 - 127.5) * dat_scl[subint, :, 3] + dat_offs[subint, :, 3]

            # Sum the corrected Stokes parameters
            summed_dynamic_spectrum = corrected_spectrum_stokes1 + corrected_spectrum_stokes2

            # Dedisperse the summed dynamic spectrum
            dedispersed_spectrum = dedisperse(summed_dynamic_spectrum, freqs, dm, sampling_time)

            # Append the dedispersed dynamic spectrum to the list
            all_summed_spectra.append(dedispersed_spectrum)

        # Concatenate all subintegrations along the time axis
        stitched_summed_spectrum = np.concatenate(all_summed_spectra, axis=0)

        # Normalize the dynamic spectrum
        stitched_summed_spectrum = (stitched_summed_spectrum - np.mean(stitched_summed_spectrum, axis=0)) / np.std(stitched_summed_spectrum, axis=0)

    # Plot the stitched and dedispersed dynamic spectrum
    time_bins = np.arange(stitched_summed_spectrum.shape[0])  # Total time bins after stitching
    frequency_channels = np.arange(stitched_summed_spectrum.shape[1])  # Frequency channels

    plt.figure(figsize=(10, 8))
    plt.imshow(stitched_summed_spectrum.T, aspect='auto', origin='lower', cmap='plasma',
               extent=[time_bins.min(), time_bins.max(), frequency_channels.min(), frequency_channels.max()],
               vmin=np.nanpercentile(stitched_summed_spectrum.T, 5),
               vmax=np.nanpercentile(stitched_summed_spectrum.T, 95))
    plt.colorbar(label='Intensity')
    plt.xlabel('Time Bins')
    plt.ylabel('Frequency Channels')
    plt.title(f'Dedispersed Dynamic Spectrum (DM = {dm})')
    plt.tight_layout()

    if output_png:
        plt.savefig(output_png)
        plt.close()
    else:
        plt.show()


def dedisperse(dynamic_spectrum, freqs, dm, sampling_time):
    """
    Apply dedispersion to the dynamic spectrum.

    Parameters:
    - dynamic_spectrum: 2D array (time x frequency) of the dynamic spectrum.
    - freqs: 1D array of frequency channels (in MHz).
    - dm: float, dispersion measure (DM) in pc cm^-3.
    - sampling_time: float, the time resolution of the data in seconds.

    Returns:
    - dedispersed_spectrum: 2D array of the dedispersed dynamic spectrum.
    """
    n_time_bins, n_freq_channels = dynamic_spectrum.shape
    dedispersed_spectrum = np.zeros_like(dynamic_spectrum)

    # Reference frequency for dedispersion (highest frequency in the band)
    ref_freq = freqs.max()

    # Compute the time delays for each frequency channel based on DM
    delays = 4.153e3 * 57.14 * (1.0 / freqs**2 - 1.0 / ref_freq**2)  # Delays in ms
    delays_in_bins = np.round(delays * 1e-3 / sampling_time).astype(int)  # Convert delays to number of time bins

    # Dedisperse by shifting each frequency channel according to its delay
    for f in range(n_freq_channels):
        shift = delays_in_bins[f]
        if shift > 0:
            dedispersed_spectrum[shift:, f] = dynamic_spectrum[:-shift, f]
        else:
            dedispersed_spectrum[:, f] = dynamic_spectrum[:, f]

    return dedispersed_spectrum


# Example usage
if __name__ == "__main__":
    plot_dedispersed_dynamic_spectrum('path_to_your_pulsar_fits_file.fits', output_png='output_dedispersed_spectrum.png')

