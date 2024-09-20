import argparse
from astropy.io import fits
import math

def fits_chunk(infname, outfname, toa):
    """
    Extract subints from a single FITS file and store them as a different FITS file.
    
    Input:
    infname: Name of input FITS file
    outfname: Name of output FITS file
    toa: Time of Arrival (TOA) used to determine subint range to extract
    
    The extracted range is based on the given TOA and the subint time (0.262144).
    """
    # Open file
    fits_file = fits.open(infname, memmap=True)

    # Read HDU 0 and 1 and their headers
    fits_hdu = fits_file[0]
    fits_hdr = fits_hdu.header

    subint_hdu = fits_file[2]
    subint_hdr = subint_hdu.header

    # Copy subint header
    new_subint_hdr = subint_hdr
    tsubint = 0.262144  # Subint time interval (seconds)
    isubmin = math.floor(toa / tsubint)
    isubmax = math.ceil(toa / tsubint + 1 / tsubint)

    print(f"Extracting subints from {isubmin} to {isubmax}")

    # Adjust NSUBOFFS in the new header
    new_subint_hdr['NSUBOFFS'] += isubmin

    # Create a new primary HDU and binary table HDU
    new_fits_hdu = fits.PrimaryHDU(data=None, header=fits_hdr)
    new_subint_hdu = fits.BinTableHDU(data=subint_hdu.data[isubmin:isubmax], header=new_subint_hdr)
    new_fits_file = fits.HDUList([new_fits_hdu, new_subint_hdu])

    # Write to new file
    new_fits_file.writeto(outfname, overwrite=True)

    # Close files
    fits_file.close()
    new_fits_file.close()

    print(f"Output written to {outfname}")
    return

# Main function to parse command line arguments
def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Extract subints from a FITS file based on Time of Arrival (TOA).")
    
    # Add arguments
    parser.add_argument('input_fits', help="Path to the input FITS file")
    parser.add_argument('output_fits', help="Path to the output FITS file")
    parser.add_argument('toa', type=float, help="Time of Arrival (TOA) to determine the subint range")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with the parsed arguments
    fits_chunk(args.input_fits, args.output_fits, args.toa)

if __name__ == "__main__":
    main()

