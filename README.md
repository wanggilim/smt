# SOFIA Dossier Generator

This package contains a series of convenience functions for accessing SOFIA DCS from the command line.  The primary purpose is to retrieve flight mission files (.mis or .misxml) from DCS, retrieve the relevant AOR files, and generate a LaTeX summary of the flight (or flight series).

Currently, the supported instruments are HAWC+, FORCAST, and FIFI-LS (partially).

## Usage

Run the dossier program as:

    sof-dossier -h

  This will print the help menu:
  
```
usage: sof-dossier [-h] [-v VERSION] [-o O] [-r] [-c] [--faor] [--sct] [--pos]
                   [-fpi] [-sio] [-preserve] [-local LOCAL] [-cfg CFG]
                   [-alias ALIAS ALIAS] [-s] [-spdf] [-se] [-reg] [-nowrite]
                   [-sf] [-irs IRSURVEY] [-mcfg MCFG] [-t TEMPLATE] [-nf]
                   flightid

Generate dossiers from flight ID

positional arguments:
  flightid              Flight ID. Can be a single flight (201809_HA_KEANU),
                        or a flight series (e.g. 201809_HA).

optional arguments:
  -h, --help            show this help message and exit
  -v VERSION, --version VERSION
                        Specify Rev number/letter
  -o O                  Output directory (default=flight series)
  -r, --refresh-cache   Force update from DCS
  -c, --compile         Compile tex to pdf
  --faor                Populate values from FAORs (FORCAST only)
  --sct                 Populate values from SCTs (FIFI-LS only)
  --pos                 Download .pos files
  -fpi, --show-fpi      Show FPI overlay
  -sio, --ioinstructions
                        Override ObsBlk comments with instrument operator
                        instructions (HAWC+ only)
  -preserve, --preserve-comments
                        Preserve comments from existing .tex file, if it
                        exists
  -local LOCAL          Specify local directory for .mis files, else query DCS
  -cfg CFG              Specify .cfg file for additional options
  -alias ALIAS ALIAS    Alias obsblocks (e.g. blk1 is blk2)
  -s, --save-aors       Save AORs to mission folder
  -spdf, --save-pdf     Save AORs and proposal pdfs to mission folder
  -se, --save-aors-exit
                        Save AORs to mission folder and exit without
                        generating dossier
  -reg, --save-regs     Save region files to mission folder
  -nowrite, --no-write  Do not [over]write .tex file.
  -sf, --save-fits      Save .fits files of fields.
  -irs IRSURVEY, --ir-survey IRSURVEY
                        Save .fits files of specified IR survey. Requires
                        --save-fits options.
  -mcfg MCFG            Model config file (default=/library/path/DBmodels.cfg)
  -t TEMPLATE, --template TEMPLATE
                        Template file (default=dossier/template.tex)
  -nf, --no-figure      Do not make figures
```


The most basic mode is the following:

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

   Finally, there are a few options that can be overridden with a config file.  An example .cfg file is attached, and can be imported as follows:

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

  The .cfg options include:

- width:  Width of the figure in degrees (default=0.4).
- height:  Height of the figure in degrees (default=0.4).
- recenter:  Specify a new center for the figure, rather than the first position in the ObsBlk.  The format is any string accepted by astropy SkyCoord (e.g. 16:20:53.50 -35:47:00.0).
- nofigure:  If 'true' or 'yes', do not display the figure for this leg/ObsBlk/planID.
- roll:  Override rotation of field angle to specified value (in degrees east of north).
- survey:  Make images from this survey (default=DSS2 Red).  This can be any survey available to astroquery.skyview.  See https://astroquery.readthedocs.io/en/latest/skyview/skyview.html for a list.
- vmin:  Override aplpy's default vmin for the image.
- vmax:  Override aplpy's default vmax for the image.
- invert: If 'false' or 'no', do not invert the colorscale (default is 'true': black on white).

