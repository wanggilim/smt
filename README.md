# SOFIA Dossier Generator

This package contains a series of convenience functions for accessing SOFIA DCS from the command line.  The primary purpose is to retrieve flight mission files (.mis or .misxml) from DCS, retrieve the relevant AOR files, and generate a LaTeX summary of the flight (or flight series).

Currently, the supported instruments are HAWC+, FORCAST, and FIFI-LS (partially).

## Usage

Run the dossier program as:

    sof-dossier -h

  This will print a help menu.  The most basic mode is the following:

    sof-dossier 201909_HA_FARUQ

   This will query DCS (asking for your DCS credentials) to find the matching flight id, download the necessary .misxml and .aor files, and generate the dossier.  It maintains an offline cache, so subsequent runs will not query DCS.  To force an update of all files related to a flight id, you can run it as:

    sof-dossier 201909_HA_FARUQ -r

  This will refresh the cache with the latest .mis and .aor files.

Running these commands will generate a folder structure in the current directory as 201909_HA/FARUQ.  Inside is the .tex file and figure pdf files.

   If you want more control over versions, it is sensible to add a version letter or number like so:

    sof-dossier 201909_HA_FARUQ -v A

  which will instead create 201909_HA_vA/FARUQ.

You can also do an entire flight series in one shot, but leaving off the specific flight name:

    sof-dossier 201909_HA -v A

This will grab all the flights associated with that flight series.

  Some other useful options include the '-nf' flag, which skips making figures (a slow step).  Also the '--compile' switch, which should attempt to automatically compile the .tex file after writing it (this requires that pdflatex is available in your PATH).  And finally the '-s' command, which saves all the .aor and proposal .pdf files to their mission folders.

  Now, there is a special case for INIT flight series that haven't made it into DCS yet.  Make sure all the init folders are grouped under one folder (an example is attached), and point the program to them with the '-local' option:

    sof-dossier 202001_HA -v A -local /path/to/init/

  This will attempt to find the .misxml (or legacy .mis) files in the local directory specified, and then query DCS for the .aor files.

   Finally, there are a few options that can be overridden with a config file.  An example is attached, and can be run as:

    sof-dossier 201909_HA -v A -cfg 201909_HA.cfg

   The section headers (in brackets) specify what the options refer to.  The section can be of an ObsBlk, planID, or aorID, applying those options to all legs in all flights with that ID. Alternatively, the section header can be specific to a certain flight, or even a certain leg.  The program will search first for the most specific criteria in reverse order as below.

- 07_0225
- OB_07_0225_01
- 07_0225_1

- GAVIN
- GAVIN_07_0225
- GAVIN_OB_07_0225_01
- GAVIN_07_0225_1

- Leg7_GAVIN
- Leg7_GAVIN_07_0225
- Leg7_GAVIN_OB_07_0225_01
- Leg7_GAVIN_07_0225_1

  The options include:
```
width:  Width of the figure in degrees (default=0.4).
height:  Height of the figure in degrees (default=0.4).
recenter:  Specify a new center for the figure, rather than the first position in the ObsBlk.  The format is any string accepted by astropy SkyCoord (e.g. 16:20:53.50 -35:47:00.0).
nofigure:  If 'true' or 'yes', do not display the figure for this leg/ObsBlk/planID.
roll:  Override rotation of field angle to specified value (in degrees east of north).
survey:  Make images from this survey (default=DSS2 Red).  This can be any survey available to astroquery.skyview.  See https://astroquery.readthedocs.io/en/latest/skyview/skyview.html for a list.
vmin:  Override aplpy's default vmin for the image.
vmax:  Override aplpy's default vmax for the image.
invert: If 'false' or 'no', do not invert the colorscale (default is 'true': black on white).
```
