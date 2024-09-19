import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
#This is a comment
#Another comment

def plot_dynamic_spectrum(fits_file, output_png=None):
    """
    Plots the stitched summed dynamic spectrum (Stokes 1 + Stokes 2) for a given FITS file,
    with DAT_SCL and DAT_OFFS applied.

    Parameters:
    - fits_file: str, path to the input FITS file.
    - output_png: str, path to save the output plot as a PNG file (optional). If None, the plot is shown interactively.
    """
    with fits.open(fits_file) as hdulist:
        # Access the binary table containing the 'DAT_SCL' and 'DAT_OFFS' columns
        data_table = hdulist[1].data  # Typically, the binary table is in the second HDU (index 1)

        # Extract the 'DAT_SCL' and 'DAT_OFFS' columns
        dat_scl = data_table['DAT_SCL']  # 'DAT_SCL' is the column name
        dat_scl = dat_scl.reshape((5, 1280,4))  # Reshape to (subint, stokes, frequency channels)
        dat_offs = data_table['DAT_OFFS']  # 'DAT_OFFS' is the column name
        dat_offs = dat_offs.reshape((5, 1280, 4))  # Reshape to (subint, stokes, frequency channels)
       # zero_offs = data_table['ZERO_OFF']
        column_names = hdulist[1].columns.names

    # Print all the column names
        print("Column names in the binary table:", column_names)
       
        print('The shape of scales is', dat_scl.shape)
        print('The shape of offsets is', dat_offs.shape)
        #print('The shape of zero_offset is', zero_offs.shape)

        # Initialize an empty list to hold the corrected dynamic spectra for each subintegration
        all_summed_spectra = []

        # Loop through each subintegration (assuming there are 5 subintegrations)
        for subint in range(5):
            # Extract the dynamic spectrum for the first two Stokes parameters
            dynamic_spectrum_stokes1 = hdulist[1].data['DATA'][subint, :, 3, :]  # First Stokes parameter (time, freq) AA
            dynamic_spectrum_stokes2 = hdulist[1].data['DATA'][subint, :, 3, :]  # Second Stokes parameter (time, freq) BB

            # Apply DAT_SCL and DAT_OFFS to both Stokes parameters
            corrected_spectrum_stokes1 = (dynamic_spectrum_stokes1-127.5) * dat_scl[subint,:,3] + dat_offs[subint,:,3]
            corrected_spectrum_stokes2 = (dynamic_spectrum_stokes2-127.5)* dat_scl[subint,:,3] + dat_offs[subint,:,3]

            # Sum the corrected Stokes parameters
            summed_dynamic_spectrum = corrected_spectrum_stokes1 + corrected_spectrum_stokes2

            # Append the summed dynamic spectrum to the list
            all_summed_spectra.append(summed_dynamic_spectrum)

        # Concatenate all subintegrations along the time axis
        stitched_summed_spectrum = np.concatenate(all_summed_spectra, axis=0)

        # Normalize the dynamic spectrum
        stitched_summed_spectrum = (stitched_summed_spectrum - np.mean(stitched_summed_spectrum, axis=0)) / np.std(stitched_summed_spectrum, axis=0)

    # Plot the stitched summed dynamic spectrum
    print(stitched_summed_spectrum.shape)
    time_bins = np.arange(stitched_summed_spectrum.shape[0])  # Total time bins after stitching
    frequency_channels = np.arange(stitched_summed_spectrum.shape[1])  # Frequency channels

    plt.figure(figsize=(10, 8))
    plt.imshow(stitched_summed_spectrum.T, aspect='auto', origin='lower', cmap='plasma',
               extent=[time_bins.min(), time_bins.max(), frequency_channels.min(), frequency_channels.max()],vmin=np.nanpercentile(stitched_summed_spectrum.T, 5), vmax=np.nanpercentile(stitched_summed_spectrum.T, 95))
    plt.colorbar(label='Intensity')
    plt.xlabel('Time Bins')
    plt.ylabel('Frequency Channels')
    plt.title('Stitched Summed Dynamic Spectrum (Stokes 1 + Stokes 2)')
    plt.tight_layout()

    if output_png:
        plt.savefig(output_png)
        plt.close()
    else:
        plt.show()


# Example usage
if __name__ == "__main__":
    plot_dynamic_spectrum('path_to_your_pulsar_fits_file.fits', output_png='output_spectrum.png')

