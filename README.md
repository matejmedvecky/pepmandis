# PepMANDIS

PepMANDIS is an automated pipeline that interrogates UniProt or user defined protein database
and computes several protein/peptide properties and associated statistics to deduce a small
list of the most representative, process-specific and MS-amenable peptides for a microbial
enzymatic activity of interest. It is written in Python3, and was tested on several Linux
distributions (Fedora, Ubuntu and CentOS) and macOS (High Sierra).

PLEASE CITE our paper and/or URL to this GitHub repo (https://github.com/matejmedvecky/pepmandis)
when using PepMANDIS in your work.

# Installation

Following list of dependencies need to be installed prior to start using pepMANDIS:

- biopython
- matplotlib
- multiprocess
- numpy
- pyopenms
- selenium

They can be installed via pip (package installer for Python). E.g. `pip install biopython`.
If, after installation of pyopenms, ImportError message (Library not loaded) is risen,
openssl need to be installed as well.
Example of installing openssl via homebrew: `brew install openssl`.

PepMANDIS pipeline uses either Safari or Chrome web browser for performing CONSeQuence analysis
(prediction of peptide LC-MS/MS detectability), therefore, one of them need to be installed
on user's computer. Please note that Chrome is the preferred one since Safari cannot be ran
in the background.\
Note for Safari users: To enable to control Safari via webdriver, 'Allow Remote Automation' option
                       in Safari's Develop menu must be enabled.\
Note for Chrome users: User system's ChromeDriver must be compatible with the system's Chrome browser
                       version. ChromeDriver can be downloaded from
                       https://sites.google.com/a/chromium.org/chromedriver/downloads.

Newer versions of macOS won't run PeptideSieve since LD_LIBRARY_PATH and DYLD_LIBRARY_PATH cannot be
loaded due to 'System integrity protection'. User must turn it off in order to enable dyld library
to be loaded.

Please note that pyOpenMS is published under 3-clause BSD licence (https://opensource.org/licenses/BSD-3-Clause).

# Quick start

`pepMANDIS.py [options] -m "desired_molecule_name"`\
e.g. `pepMANDIS.py -m "catechol-1,2-dioxygenase"`

If path to pepMANDIS.py is not in your environment PATH, execute it as follows:\
`/path/to/pepMANDIS.py [options] -m "desired_molecule_name"`

If the script is in your current working directory:\
`./pepMANDIS.py [options] -m "desired_molecule_name"`

Help:\
`pepMANDIS -h`

# Licence

[GPL Version 3](https://github.com/matejmedvecky/pepmandis/blob/master/LICENSE)

# Author

[Matej Medvecky](https://github.com/matejmedvecky)

# Credits

Matej Medvecky\
Manolis Mandalakis
