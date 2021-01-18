# PepMANDIS

PepMANDIS is an automated pipeline that interrogates UniProt or user defined protein database
and computes several protein/peptide properties and associated statistics to deduce a small
list of the most representative, process-specific and MS-amenable peptides for a microbial
enzymatic activity of interest. It is written in Python, and was tested on several Linux
distributions (Fedora, Ubuntu and CentOS), macOS (High Sierra) and Windows (10).

PLEASE CITE our paper (currently under revision) and/or URL to this GitHub repo (https://github.com/matejmedvecky/pepmandis)
when using PepMANDIS in your work.

# Installation

Following list of dependencies need to be installed prior to start using pepMANDIS:

- python 3 (3.7.6 High Sierra, 3.8.5 Windows 10)
- biopython (1.76 High Sierra, 1.78 Windows 10)
- matplotlib (3.1.3 High Sierra, 3.3.3 Windows 10)
- multiprocess (0.70.9 High Sierra, 0.70.11 Windows 10)
- numpy (1.18.1 High Sierra, 1.19.5 Windows 10)
- pyopenms (2.4.0 High Sierra, 2.6.0 Windows 10)
- selenium (3.141.0 High Sierra, Windows 10)



- BLAST+ (2.10.0+ macOS) & BLAST non-redundant protein database (HIGHLY RECOMMENDED)
- PeptideSieve (0.51)

They can be installed via pip (package installer for Python) in Terminal. E.g. `pip install biopython`.
If, after installation of pyopenms, ImportError message (Library not loaded) is risen,
openssl need to be installed as well.
Example of installing openssl via homebrew in Terminal: `brew install openssl`.

PepMANDIS pipeline uses either Safari or Chrome web browser for performing CONSeQuence analysis
(prediction of peptide LC-MS/MS detectability), therefore, one of them need to be installed
on user's computer. Please note that Chrome is the preferred one since Safari cannot be ran
in the background.\
Note for Safari users: To enable to control Safari via webdriver, 'Allow Remote Automation' option
                       in Safari's Develop menu must be enabled.\
Note for Chrome users: User system's ChromeDriver must be compatible with the system's Chrome browser
                       version. ChromeDriver can be downloaded from
                       https://sites.google.com/a/chromium.org/chromedriver/downloads.

Non-redundant BLAST database can be downloaded via following FTP:  https://ftp.ncbi.nlm.nih.gov/blast/db/. \
Please note that path to BLAST NR protein DB need to be specified in a 'defaults.cfg' file. 

Newer versions of macOS won't run PeptideSieve since LD_LIBRARY_PATH and DYLD_LIBRARY_PATH cannot be
loaded due to 'System integrity protection'. User must turn it off in order to enable dyld library
to be loaded.

Please note that pyOpenMS is published under 3-clause BSD licence (https://opensource.org/licenses/BSD-3-Clause).

# Quick start

`pepMANDIS.py [options] -m "desired_molecule_name"`\
e.g. `pepMANDIS.py -m "catechol-1,2-dioxygenase"`\
or  `python pepMANDIS.py -m "catechol-1,2-dioxygenase"`

Help:\
`pepMANDIS -h`

## Linux and macOS

Open the terminal (press Ctrl+Alt+T in Linux or press Command+Space, type 'terminal' and press Enter in macOS).

If path to pepMANDIS.py is not in your environment PATH, execute it as follows:\
`/path/to/pepMANDIS.py [options] -m "desired_molecule_name"`

If the script is in your current working directory, execute it by following way:\
`./pepMANDIS.py [options] -m "desired_molecule_name"`

## Windows

Open the terminal (press Windows+R, type 'cmd' and press Enter).

If path to pepMANDIS.py is not in your environment PATH, execute it as follows:\
`C:\path\to\pepMANDIS.py [options] -m "desired_molecule_name"`\
or `python C:\path\to\pepMANDIS.py -m "catechol-1,2-dioxygenase"`

If the script is in your current working directory, execute it by following way:\
`pepMANDIS.py [options] -m "desired_molecule_name"`\
or `python pepMANDIS.py -m "catechol-1,2-dioxygenase"`

# Licence

[GPL Version 3](https://github.com/matejmedvecky/pepmandis/blob/master/LICENSE)

# Author

[Matej Medvecky](https://github.com/matejmedvecky)

# Credits

Matej Medvecky\
Manolis Mandalakis
