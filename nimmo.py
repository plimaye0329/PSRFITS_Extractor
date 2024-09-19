from astropy.io import fits
import math

def fits_chunk(infname,outfname,toa):
    """Extract subints from a single FITS file and store them as a different FITS file

    Input:
       infname: Name of input FITS file
       outfname: Name of output FITS file
       isubmin,isubmax: subint range to extract and store

    """
    # Open file
    fits_file=fits.open(infname,memmap=True)

    # Read HDU 0 and 1 and their headers
    fits_hdu=fits_file[0]
    fits_hdr=fits_hdu.header
    #hdr = fits_file[0].header
    #scl = hdr['DAT_SCL']
    #offs = hdr['DAT_OFFS']
    #print('Data')
    subint_hdu=fits_file[2]
    subint_hdr=subint_hdu.header

   

    # Copy subint header
    new_subint_hdr=subint_hdr
    tsubint = 0.262144
    isubmin = math.floor(toa/tsubint)
    isubmax = math.ceil(toa/tsubint + 1/tsubint)

    print(isubmin)
    print(isubmax)

    

    # Adjust NSUBOFFS
    new_subint_hdr['NSUBOFFS']+=isubmin

    # Create a new primary HDU and binary table HDU
    new_fits_hdu=fits.PrimaryHDU(data=None,header=fits_hdr)
    new_subint_hdu=fits.BinTableHDU(data=subint_hdu.data[isubmin:isubmax],header=new_subint_hdr)
    new_fits_file=fits.HDUList([new_fits_hdu,new_subint_hdu])

    # Write to new file
    new_fits_file.writeto(outfname,overwrite=True)

    # Close files
    fits_file.close()
    new_fits_file.close()

    return 
