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

## Mandatory dependencies

- python 3 (3.7.6 High Sierra, 3.8.5 Windows 10)
- biopython (1.76 High Sierra, 1.78 Windows 10)
- matplotlib (3.1.3 High Sierra, 3.3.3 Windows 10)
- multiprocess (0.70.9 High Sierra, 0.70.11 Windows 10)
- numpy (1.18.1 High Sierra, 1.19.5 Windows 10)
- pyopenms (2.4.0 High Sierra, 2.6.0 Windows 10)
- selenium (3.141.0 High Sierra, Windows 10)

They can be installed via pip (package installer for Python) in Terminal. E.g. `pip install biopython`.
If, after installation of pyopenms, ImportError message (Library not loaded) is risen,
openssl needs to be installed as well.
Example of installing openssl via homebrew in Terminal: `brew install openssl`.

Terminal can be opened by following ways:\
Linux: press Ctrl+Alt+T\
macOS: press Command+Space, type `terminal` and press Enter\
Windows: press Windows+R, type `cmd` and press Enter

PepMANDIS pipeline uses either Safari or Chrome web browser for performing CONSeQuence analysis
(prediction of peptide LC-MS/MS detectability), therefore, one of them need to be installed
on user's computer. Please note that Chrome is the preferred one since Safari cannot be ran
in the background.\
Note for Safari users: To enable to control Safari via webdriver, 'Allow Remote Automation' option
                       in Safari's Develop menu must be enabled.\
Note for Chrome users: User system's ChromeDriver must be compatible with the system's Chrome browser
                       version. ChromeDriver can be downloaded from
                       https://sites.google.com/a/chromium.org/chromedriver/downloads.

Please note that pyOpenMS is published under 3-clause BSD licence (https://opensource.org/licenses/BSD-3-Clause).

## Optional dependencies

- BLAST+ (2.10.0+ High Sierra) & BLAST non-redundant protein database (HIGHLY RECOMMENDED)
- PeptideSieve (0.51)

Non-redundant BLAST database can be downloaded via following FTP:  https://ftp.ncbi.nlm.nih.gov/blast/db/. \
Please note that path to BLAST NR protein DB need to be specified in a 'defaults.cfg' file. 

Newer versions of macOS won't run PeptideSieve since LD_LIBRARY_PATH and DYLD_LIBRARY_PATH cannot be
loaded due to 'System integrity protection'. User must turn it off in order to enable dyld library
to be loaded.

## PepMANDIS source code

PepMANDIS repository with the souce code can be obtained via git:\
`git clone https://github.com/matejmedvecky/pepmandis.git`\
Scripts are located in bin directory:\
They can by accessed on Linux/macOS:\
`cd pepmandis/bin`

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
or `python C:\path\to\pepMANDIS.py -m "desired_molecule_name"`

If the script is in your current working directory, execute it by following way:\
`pepMANDIS.py [options] -m "desired_molecule_name"`\
or `python pepMANDIS.py -m "desired_molecule_name"`

# Options

- `-m`, `--molecule`: Specify protein name bounded by quotes (e.g. "catechol 1,2-dioxygenase") \
  that is used to programmatically query UniProt database and retrieve all relevant \
  entries, and to generate a list of potential protein synonyms reported in BLAST \
  non-redundant protein DB. REQUIRED.

- `--no-usearch`: Use this option if input sequences are provided in a FASTA AA \
  file (argument `-i`). Gathering of UniProt data will be avoided.

- `--extra-input`: Use this option if an additional dataset with protein sequences in FASTA AA \
  format is also provided (argument `-e`). Such proteins will be compared to \
  entries gathered from UniProt DB or sequences provided by '-i' argument \
  and exact matching peptides only will be used in the downstream analysis. \
  (See manual or manuscript for more detailed info).

- `--custom-url`: UniProt search URL is provided as an argument; program will \
  extract input sequences using provided URL. `--utaxonomy` arg is omitted \
  if this argument is specified.

- `--no-bsearch`: Do not perform BLASTP search. Use this option also if using online blasting and \
  'blastp.xml' file has already been generated in previous run in order \
  to save considerable amout of time.

- `--offline-blastp`: Perform offline blasting (RECOMMENDED). Provide path to nr database (databasePath) in a \
  config file 'defaults.cfg'.

- `--run-peptidesieve`: By default, CONSeQuence tool is only run in order to predict detectability \
  of peptides by MS. This option enables PeptideSieve to be executed as well.

- `-C`, `--config`: Specify path to the configuration file (e.g. /Users/me/Documents/defaults.cfg). \
  [default: ./defaults.cfg (current working directory)].

- `-l`, `--length`: Coefficient of allowed length variation of entries. Allowed \
				'lengths interval = (median length value) +- (coefficient of l \
				'variation)\*(median length value). Entries with higher length \
				'difference are removed. [default: 0.25]

- `--no-putatives`: Use this option if putative molecules from UniProt search should \
  be removed.

- `-u`, `--utaxonomy`: Restrict protein search in UniProt database for organisms \
  specified by taxonomy keywords bounded by quotes: for AND use \
  character '&', for OR use character '|', for NOT use character \
  '-'. (E.g. "Actinobacteria|Proteobacteria-Pseudomonadales" which \
  means Actinobacteria OR Proteobacteria without (NOT) \
  Pseudomonadales). [default: "Bacteria"]

- `-c`, `--cleavages`: Number of missed cleavages for the trypsin digestion of proteins. \
  [default: 0]

- `-x`, `--mins`: Smallest size [in AA] of digested peptide to be kept in a \
  dataset. [default: 8]

- `-y`, `--maxs`: Largest size [in AA] of digested peptide to be kept in a dataset. \
  [default: 25]

- `-b`, `--btaxonomy`: Limit BLASTP search against proteins belonging to organisms \
  specified by taxonomy keywords bounded by quotes: for AND use \
  character '&', for OR use character '|', for NOT use character \
  '-'. (E.g. "Actinobacteria|Proteobacteria-Pseudomonadales" which \
  means Actinobacteria OR Proteobacteria without (NOT) \
  Pseudomonadales). [default: "Bacteria"]

- `-k`,`--ktaxonomy`: Specified genera/species (genera is represented by just one \
  word; species is represented by just two words separated by space) \
  will be represented by 3 peptides (per genera/species) in a final \
  list of selected peptides. Bound specified taxonomy by quotes; for \
  AND use character '&' (E.g. "Pseudomonas&Pseudomonas putida" will \
  output 3 peptides representing genera Pseudomonas in file \
  'Selected_peptides.genera' as well as 3 peptides representing \
  species Pseudomonas putida in file 'Selected_peptides.species'). \
  [default: None]

- `-s`, `--specificity`: Peptide specificity score threshold. [default: 90.0]

- `-a`, `--avoid-spec-filt-list`: Provide a list of peptides for that specificity scores are omitted. \
  I.e. specificity score-based filtering step will not apply for such peptides. \
  Bound specified peptides by quotes and separate individual entries by comma \
  (E.g. "ENPPVLPK,SGLFTSEELPR").

- `-cd`, `--cdetectability`: Minimum number of CONSeQence ML algorithms (out of 4) considering a peptide \
  as detectable by MS to pass detectability filter. Please note that maximum \
  value is 4. [default: 1]

- `-pd`, `--pdetectability`: PeptideSieve detectability score threshold. [default: 0.6]

- `-f`, `--chemmod-filter`: Minimum number of possible chemical modifications of a peptide \
  to be omitted from final selection lists. [default: 4]

- `-n`, `--peptide-calc-count`: Number of peptides with the highest coverage in terms of microbial species \
  to calculate specificity and detectability scores for. [default: 400]

- `-N`, `--peptide-out-count`: Number of peptides to be output to 'Selected_peptides.\*' files. [default: 50]

- `-t`, `--threads`: Number of threads (processes) to be used in parallel blasting. \
  Works for online as well as offline blasting. [default: 4]

- `-T`, `--timeout`: Maximum time (in seconds) to wait for online BLASTP results. \
  [default: 1500]

- `-z`, `--min-synname-length`: Minimum number of words blastp synonymous name consist of to be \
  accepted as synonymous name. [default: 1]

- `-U`, `--in-url`: Specify UniProt URL in quotes. Applicable only if `--custom-url` \
  argument is specified.

- `-i`, `--in-aa-file`: Specify path to dataset with protein sequences in FASTA AA format (e.g. \
  infile.faa). Applicable only if `--no-usearch` argument is specified. Following  \
  header formats are allowed: >PROTEIN_ID or >PROTEIN_ID~COV=<value> if \
  coverage (expected copy number) of entries is also provided. \
  Value can be either integer or float larger than 1.0.

- `-e`, `--in-extra-aa-file`: Specify path to an extra input fasta AA file (e.g. infile_extra.faa). \
  Applicable only if `--extra-input` argument is specified. Following \
  header formats are allowed: >PROTEIN_ID or >PROTEIN_ID~COV=<value> in \
  case of coverage (expected copy number) of entries is also provided. \
  Value can be either integer or float larger than 1.0.

- `-I`, `--in-blastp-file`: Specify path to 'blastp.xml' file generated in previous program \
  run. Applicable only if `--no-bsearch` argument is specified.

- `--stats-only`: Print general info, protein statistics and figures, and exit.

# Licence

[GPL Version 3](https://github.com/matejmedvecky/pepmandis/blob/master/LICENSE)

# Author

[Matej Medvecky](https://github.com/matejmedvecky)

# Credits

Matej Medvecky\
Manolis Mandalakis

# Citation

Medvecky M, Mandalakis M. PepMANDIS: A peptide selection tool for designing function-based targeted proteomic assays in mixed microbial communities. *Currently under revision.*
