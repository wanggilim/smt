# SOFIA Mission Toolbox

This package contains a series of convenience functions for accessing SOFIA DCS from the command line.  The primary purpose is to retrieve flight mission files (.mis or .misxml) from DCS, retrieve the relevant AOR files, and generate a LaTeX summary of the flight (or flight series).

Currently, the supported instruments are HAWC+, FORCAST, and FIFI-LS (partially).

The architecture of the following programs was designed initially as a means of automating common USPOT/DCS interactions.  Most of the programs use a shared database locally to cache any and all files required from DCS.  Should the need arise, dossiers can be generated on-the-fly with no internet connectivity (i.e. on a SOFIA flight), as long as the relevant DCS data has been previously cached.  Much of this caching/database interaction is handled automatically and invisibly from the user.

## Installation
To install the dossier generator and assorted other programs, clone this repository to your machine and run the following command:

	python setup.py install
	
The necessary dependencies should download and install automatically, but if not, read whatever error messages appear and try to track down the required packages manually, either through `conda install` or `pip install`.

If you have previously installed this toolbox, run `sof-refresh` immediately following the install.

## Programs

The install procedure should add all of the following programs to your `PATH` variable:

- [`sof-dossier`](#markdown-header-sof-dossier):  Generate dossiers from flight/series ID.
- [`sof-planner`](#markdown-header-sof-planner):  Generate FORCAST `.faor` files.
- [`aor2faor`](#markdown-header-aor2faor):  Convert local `.aor` files to `.faor` for FORCAST.
- [`sof-itime`](#markdown-header-sof-itime):  FORCAST integration time efficiency calculator.
- [`sof-aor`](#markdown-header-sof-aor):  Utility program to download `.aor` and proposal `.pdf` files to a local directory.
- [`sof-cache`](#markdown-header-sof-cache):  Utility program for pre-caching all necessary files for `sof-dossier`, `sof-planner`, or `sof-aor`.
- [`sof-refresh`](#markdown-header-sof-refresh):  Force refresh of local database from previously-cached DCS files.

Documentation for each of these programs is described below.  All of them have a help menu that appears with the addition of the `-h` command line switch.


## sof-dossier

Generate dossiers from flight/series ID

### Usage

Run the dossier program as:

    sof-dossier -h

```bash
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

  This will refresh the cache with the latest `.misxml` and `.aor` files.

Running these commands will generate a folder structure in the current directory as 201909_HA/FARUQ.  Inside are the `.tex` file and figure `.pdf` files for the dossier.

   If you want more control over versions, it is sensible to add a version letter or number like so:

    sof-dossier 201909_HA_FARUQ -v A

  which will instead create 201909_HA_vA/FARUQ.

You can also do an entire flight series in one shot by leaving off the specific flight name:

    sof-dossier 201909_HA -v A

This will grab all the flights associated with that flight series.

  Some other useful options include the `-nf` flag, which skips making figures (a slow step).  Also the `--compile` switch, which should attempt to automatically compile the `.tex` file after writing it (this requires that pdflatex is available in your PATH).  And finally the `-s` command, which saves all the `.aor` and proposal .pdf files to their mission folders.

  Now, there is a special case for INIT flight series that haven't made it into DCS yet.  Make sure all the init folders are grouped under one folder (an example is attached), and point the program to them with the `-local` option:

    sof-dossier 202001_HA -v A -local /path/to/init/

  This will attempt to find the `.misxml` (or legacy `.mis`) files in the local directory specified, and then query DCS for the .aor files.

   Finally, there are a few options that can be overridden with a config file.  An example `.cfg` file is attached, and can be imported as follows:

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

  The `.cfg` options include:

- width:  Width of the figure in degrees (default=0.4).
- height:  Height of the figure in degrees (default=0.4).
- recenter:  Specify a new center for the figure, rather than the first position in the ObsBlk.  The format is any string accepted by astropy SkyCoord (e.g. 16:20:53.50 -35:47:00.0).
- nofigure:  If 'true' or 'yes', do not display the figure for this leg/ObsBlk/planID.
- roll:  Override rotation of field angle to specified value (in degrees east of north).
- survey:  Make images from this survey (default=DSS2 Red).  This can be any survey available to astroquery.skyview.  See https://astroquery.readthedocs.io/en/latest/skyview/skyview.html for a list.
- vmin:  Override aplpy's default vmin for the image.
- vmax:  Override aplpy's default vmax for the image.
- invert: If 'false' or 'no', do not invert the colorscale (default is 'true': black on white).

To add in additional information from `.faor` files (if created with the `sof-planner` program described below), include the `--faor` command line switch.


## sof-planner

Generate FORCAST `.faor` files.

### Usage

Run the planner program as:

    sof-planner -h

```
usage: sof-planner [-h] [-niter NITER] [--basic] [--quick]
                   [-losbuffer LOSBUFFER] [--debug] [-leg LEG]
                   [-interval INTERVAL] [--allaors] [-alias ALIAS ALIAS]
                   [-rofrate ROFRATE] [--dry-run] [-o O] [-v VERSION]
                   [-local LOCAL] [-cfg CFG] [-mcfg MCFG] [--comment-out]
                   [--fixed] [-noddwell NODDWELL] [-loops LOOPS]
                   [-repeats REPEATS] [-dithers DITHERS] [-r]
                   flightid

Plan FORCAST observations, accounting for LOS rewind cadence.

positional arguments:
  flightid              Flight ID (e.g. 201902_FO_OTTO)

optional arguments:
  -h, --help            show this help message and exit
  -niter NITER          Number of iterations per AORID (default=6)
  --basic               If specified, generate basic FAORs without
                        optimization
  --quick               If specified, generate FAORs with quick optimization.
                        May not find optimal solution.
  -losbuffer LOSBUFFER  Allowed time in sec to go over LOS rewind cadence
                        (default=30)
  --debug               If specified, disable multiprocessing
  -leg LEG              Only process this leg
  -interval INTERVAL    Split leg into intervals of this many minutes
                        (requires -leg)
  --allaors             If specified, ignore aorids listed in config
  -alias ALIAS ALIAS    Alias obsblocks (e.g. blk1 is blk2)
  -rofrate ROFRATE      Fix rof rate (in deg/min) for all legs
  --dry-run             If specified, print plan and exit
  -o O                  Output directory (default=flight series)
  -v VERSION, --version VERSION
                        Specify Rev number/letter
  -local LOCAL          Specify local directory for .mis files, else query DCS
  -cfg CFG              Specify .cfg file for additional options
  -mcfg MCFG            Model config file (default=/library/path/DBmodels.cfg)
  --comment-out, -co    If specified, leave unplanned run blocks in FAORs
                        rather than delete them
  --fixed               If specified, do not iterate, and use initial values.
  -noddwell NODDWELL    Specify noddwell for all legs
  -loops LOOPS          Specify loops for all legs
  -repeats REPEATS      Specify repeats for all legs
  -dithers DITHERS      Specify dithers for all legs
  -r, --refresh-cache   Force update from DCS
```

### Description

When run, this program downloads `.aor` and `.misxml` files from DCS (unless already stored in the local cache from a previous `sof-planner` or `sof-dossier` run), and attempts to automatically produce corresponding `.faor` files for each leg in a flight/flight series.  Once complete, the configurations in the `.faor` files can be incorporated into the dossiers by running `sof-dossier` with the `--faor` switch.

Many of the arguments for this program are the same as `sof-dossier`.  Notably, the only required argument is a flight series, or single flight plan.  Perhaps obviously, the flight ID should correspond toa FORCAST flight series (_FO).  Similarly to the dossier program, the `-local` flag can be used to point to series `.misxml` files that have been downloaded locally, rather than pulling from DCS.  Also, the `-r` flag can be used to pull fresh `.aor` and `.misxml` files from DCS and update the local cache.

The automated fitting for optimal FORCAST configurations can be overridden on the command line with the `-nodwell`, `-loops`, `-repeats`, and `-dither` arguments (and others) followed by the desired value.  Note that these options apply globally to every leg processed.  To restrict an optimization/`.faor` generation run to just a single leg (either to override the values on the command line, or just to save time), add the `-leg` switch ti select just a single leg/program in the flight.

If you make changes to `.aors` in USPOT, and you want to update both the `.faor` and dossier files, you can always re-run the above commands with `-r` to refresh the cache.  You only need to do the refresh with one of the commands (so, do it with `sof-planner`).  Again, if you only made changes to one program, you can also add `-r -leg 7` to `sof-planner` to only generate new faors for, in this case, leg  7, which will save some time and will not overwrite the files for the other legs.

Just as in `sof-dossier` each leg/project/aorID can be customized with a `.cfg` file.  An example FORCAST `.cfg` file is provided in this repo.

To generate just the basic `.faor` files with no optimization, add the `--basic` flag.  This roughly reproduces the behavior of Kaori Nishikada's program, as well as the `aor2faor` program described below.

To manually tweak a configuration, I would recommend using the `sof-itime` program described below.  This recreates the IDL `itime3` program from Jim De Buizer, and allows for estimating the total duration of noddwell, repeat, loop, and dither combinations (for NMC, C2NC2, and NXCAC modes).  When desired values are found, you can add the appropriate noddwell, repeat, etc. parameters to a `.cfg` file and re-run `sof-planner`.


## aor2faor

Generate basic FORCAST `.faor` files

### Usage

Run the aor translator program as:

    aor2faor -h

```
usage: aor2faor [-h] [-o outdir] aors [aors ...]

AOR to FAOR translator

positional arguments:
  aors        .aor file(s) to process

optional arguments:
  -h, --help  show this help message and exit
  -o outdir   Output directory (default='.')
```

### Description

This program converts `.aor` files with FORCAST AORs into basic `.faor` files.  The output corresponds to `sof-planner` with the `--basic` flag, i.e. FAORs with no optimization.


## sof-itime

FORCAST integration time efficiency calculator

### Usage

Run the calculator as:

	sof-itime -h
	
```
usage: sof-itime [-h] [-loops LOOPS] [-chopeff CHOPEFF] [--c2nc2] [--nxcac]
                 [-nodtime NODTIME] [-dithertime DITHERTIME]
                 [-filtchange FILTCHANGE] [-lost LOST]
                 noddwell repeats [rewinds] [dithers]

Integration time efficiency calculator

positional arguments:
  noddwell              Nod dwell time in sec
  repeats               Number of repeats
  rewinds               Number of repeats per rewind e.g. rewind cadence
                        (default=0; no rewinds during AOR)
  dithers               Number of dithers (default=0)

optional arguments:
  -h, --help            show this help message and exit
  -loops LOOPS          Number of loops for C2NC2 mode (default=1)
  -chopeff CHOPEFF      Chop efficiency (default=0.70)
  --c2nc2               C2NC2 mode
  --nxcac               NXCAC mode
  -nodtime NODTIME      Nod overhead in sec (default=5)
  -dithertime DITHERTIME
                        Dither overhead in sec (default=5)
  -filtchange FILTCHANGE
                        Filter change overhead in sec (default=30)
  -lost LOST            Additional overhead in sec (default=25)
```

### Description

This program allows for manually iterating through a sequence of noddwells, repeats, rewinds, etc. to achieve a desired total integration time or total AOR execution time.  This is helpful for planning observations where ROF rates/leg lengths are a concern on timing.  The `sof-planner` program uses the timing and efficiency from this calculator to estimate an optimal solution.

## sof-aor

Utility program to download `.aor` and proposal `.pdf` files to a local directory.

### Usage

Run the utility as:

	sof-aor -h

```
usage: sof-aor [-h] [-o O] [-pdf] [-r] planids [planids ...]

Download AORs locally

positional arguments:
  planids              Input planID(s)

optional arguments:
  -h, --help           show this help message and exit
  -o O                 Output directory (default=.)
  -pdf                 Download proposal pdfs as well
  -r, --refresh-cache  Force update from DCS
```

### Description

This program pulls `.aor` (and optionally the program proposal) files from DCS for a given planID.  It first attempts to copy the file locally from the cache, if previously downloaded through `sof-dossier`, `sof-planner`, or `sof-cache`, to the desired folder.  If the `.aor` file has not been downloaded previously (or if the `-r` option is specified), the program will pull the file directly from DCS and update the local cache.

This utility is a quick way of pulling AORs without using USPOT or DCS.

## sof-cache

Utility program for pre-caching all necessary files for `sof-dossier`, `sof-planner`, or `sof-aor`.

### Usage

Run the utility as:

	sof-cache -h

```
usage: sof-cache [-h] [-local LOCAL] [-alias ALIAS ALIAS] [-mcfg MCFG]
                 ids [ids ...]

Refresh cache for specified flight/series/obsplan.

positional arguments:
  ids                 ObsPlan IDs, flight plans, or series

optional arguments:
  -h, --help          show this help message and exit
  -local LOCAL        Specify local directory for .mis files, else query DCS
  -alias ALIAS ALIAS  Alias obsblocks (e.g. blk1 is blk2)
  -mcfg MCFG          Model config file (default=/library/path/DBmodels.cfg)
```

### Description

This program is effectively just the back-end/precursor to `sof-dossier` or `sof-planner`.  When a flight series, flight plan ID, or ObsPlanID(s) are given, the program pulls all relevant files from DCS to the local database.  If the same flight plan, flight series, or planIDs are requested from the other programs above, those programs should run significantly faster, as no new files will be needed from DCS.

## sof-refresh

Force refresh of local database from previously-cached DCS files

### Usage

Run this helper utility as:

	sof-refresh -h

```
usage: sof-refresh [-h] [-mcfg MCFG]

Force refresh of smt.db file from local cache

optional arguments:
  -h, --help  show this help message and exit
  -mcfg MCFG  Model config file (default=/library/path/DBmodels.cfg)
```

### Description

This program should be run whenever the SOFIA Mission Toolbox software is updated (in particular, when modifications are made to the database schema).  This may also help with other issues/error messages.

## sof-obsmaker

(BETA) Generate FIFI sct and map files.

### Usage

Run this program as:

	sof-obsmaker -h

```
usage: sof-obsmaker [-h] [-leg LEG] [-alias ALIAS ALIAS] [-o O] [-v VERSION]
                    [-local LOCAL] [-ocfg OCFG] [-mcfg MCFG] [-r]
                    flightid

Generate FIFI sct and map files.

positional arguments:
  flightid              Flight ID (e.g. 201902_FO_OTTO)

optional arguments:
  -h, --help            show this help message and exit
  -leg LEG              Only process this leg
  -alias ALIAS ALIAS    Alias obsblocks (e.g. blk1 is blk2)
  -o O                  Output directory (default=flight series)
  -v VERSION, --version VERSION
                        Specify Rev number/letter
  -local LOCAL          Specify local directory for .mis files, else query DCS
  -ocfg OCFG            Obsmaker config file
                        (default=/library/path/obsmaker.cfg)
  -mcfg MCFG            Model config file (default=/library/path/DBmodels.cfg)
  -r, --refresh-cache   Force update from DCS
```

Similar to `sof-planner`, this program produces the FIFI configuration files.  These can be incorporated into the dossiers with the `sof-dossier --sct` switch.
