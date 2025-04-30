## Extract bursts into 2second FITS snapshots:
This repo consists of scripts and examples of PSRFITS extractions.

The Ultra BroadBand receiver covers a frequency range of 1.3 - 6 GHz divided into five sub-bands. The FRB observations therefore demand a heavy data storage.
The PSRFITS Extraction script performs offline burst extraction post their detection using TransientX. FRB observations are mostly noise except for a few detections. Hence, after the extraction, one could erase the data to save diskspace. Alternatively, the original data could be saved to tape archives and could be recovered in case a periodicity is discovered in any observed FRB :-)
