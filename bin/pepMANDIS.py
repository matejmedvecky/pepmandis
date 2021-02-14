#!/usr/bin/env python

"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import re
import sys
import os
import pathlib
from pyopenms import AASequence, ProteaseDigestion # pyopenms has issues with Qt*dll files in Windows, therefore, it must be imported before other specific libraries (such as matplotlib)
import urllib.request, urllib.parse, urllib.error
import traceback
import numpy
import matplotlib.pyplot
import multiprocess
import time
import subprocess
import glob
import platform
from time import localtime, strftime
from configparser import ConfigParser
from itertools import chain
from collections import Counter, OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from selenium import webdriver

__author__ = 'Matej Medvecky'
__copyright__ = 'Copyright 2021, Matej Medvecky'
__credits__ = ['Matej Medvecky', 'Manolis Mandalakis']
__license__ = 'GPLv3'
__maintainer__ = 'Matej Medvecky'
__email__ = 'mat.medvecky@gmail.com'
__status__ = 'Production'
__version__ = '1.0'


class MetaproteomicPipline:

	def __init__(self):
		self.proteinEntries = []# each entry consist of protein id, coverage, organism taxonomy, phylum, length, aa sequence, digested peptides
		self.resultsFolder = ''
		self.phylumEntries = []
		self.peptideEntries = []# each entry consist of aa sequence, frequency, different species count, peptideSieve detectability score, CONSeQuence detectability, specificity score, possible chemical modifications count
		self.queryPeptides = []

	def parse_user_input(self):
		parser = argparse.ArgumentParser(
			description=
			'PepMANDIS is an automated pipeline that interrogates UniProt or user defined protein ' \
			'database and computes several protein/peptide properties and associated statistics ' \
			'to deduce a small list of the most representative, process-specific and MS-amenable ' \
			'peptides for a microbial enzymatic activity of interest.'
			)
		parser.add_argument('-m', '--molecule',
			help='Specify protein name bounded by quotes (e.g. "catechol 1,2-dioxygenase") ' \
				'that is used to programmatically query UniProt database and retrieve all relevant ' \
				'entries, and to generate a list of potential protein synonyms reported in BLAST ' \
				'non-redundant protein DB. REQUIRED.',
			type=str, required=True, metavar="\"MOLECULE\"",
			)
		parser.add_argument('--no-usearch',
			help='Use this option if input sequences are provided in a FASTA AA ' \
				'file (argument \'-i\'). Gathering of UniProt data will be avoided.',
			required=False, action='store_true', default=False,
			)
		parser.add_argument('--extra-input',
			help='Use this option if an additional dataset with protein sequences in FASTA AA ' \
				'format is also provided (argument \'-e\'). Such proteins will be compared to ' \
				'entries gathered from UniProt DB or sequences provided by \'-i\' argument ' \
				'and exact matching peptides only will be used in the downstream analysis. ' \
				'(See manual or manuscript for more detailed info)',
			required=False, action='store_true', default=False,
			)
		parser.add_argument('--custom-url',
			help='UniProt search URL is provided as an argument; program will ' \
				'extract input sequences using provided URL. --utaxonomy arg is omitted ' \
				'if this argument is specified.',
			required=False, action='store_true', default=False,
			)
		parser.add_argument('--no-bsearch',
			help='Do not perform BLASTP search. Use this option also if using online blasting and ' \
				'\'blastp.xml\' file has already been generated in previous run in order ' \
				'to save considerable amout of time.',
			required=False, action='store_true', default=False,
			)
		parser.add_argument('--offline-blastp',
			help='Perform offline blasting (RECOMMENDED). Provide path to nr database (databasePath) in a ' \
			'config file \'defaults.cfg\'.',
			required=False, action='store_true', default=False
			)
		parser.add_argument('--run-peptidesieve',
			help='By default, CONSeQuence tool is only run in order to predict detectability ' \
				'of peptides by MS. This option enables PeptideSieve to be executed as well.',
			required=False, action='store_true', default=False
			)
		parser.add_argument('-C', '--config',
			help='Specify path to the configuration file (e.g. /Users/me/Documents/defaults.cfg).' \
				'[default: ./defaults.cfg (current working directory)]',
			type=str, required=False, default="./defaults.cfg",
			)
		parser.add_argument('-l', '--length',
			help='Coefficient of allowed length variation of entries. Allowed ' \
				'lengths interval = (median length value) +- (coefficient of l ' \
				'variation)*(median length value). Entries with higher length ' \
				'difference are removed. [default: 0.25]',
			type=float, required=False, default=0.25, metavar="L_VAR_COEF",
			)
		parser.add_argument('--no-putatives',
			help='Use this option if putative molecules from UniProt search should ' \
			'be removed.',
			required=False, action='store_true', default=False,
			)
		parser.add_argument('-u', '--utaxonomy',
			help='Restrict protein search in UniProt database for organisms ' \
				'specified by taxonomy keywords bounded by quotes: for AND use ' \
				'character \'&\', for OR use character \'|\', for NOT use character ' \
				'\'-\'. (E.g. "Actinobacteria|Proteobacteria-Pseudomonadales" which ' \
				'means Actinobacteria OR Proteobacteria without (NOT) ' \
				'Pseudomonadales). [default: "Bacteria"]',
			type=str, required=False, default="Bacteria",
			metavar="\"UNIPROT_TAXONOMY\"",
			)
		parser.add_argument('-c', '--cleavages',
			help='Number of missed cleavages for the trypsin digestion of proteins. ' \
				'[default: 0]',
			type=int, required=False, default=0,
			)
		parser.add_argument('-x', '--mins',
			help='Smallest size [in AA] of digested peptide to be kept in a ' \
				'dataset. [default: 8]',
			type=int, required=False, default=8,
			)
		parser.add_argument('-y', '--maxs',
			help='Largest size [in AA] of digested peptide to be kept in a dataset. ' \
				'[default: 25]',
			type=int, required=False, default=25,
			)
		parser.add_argument('-b', '--btaxonomy',
			help='Limit ONLINE BLASTP search against proteins belonging to organisms ' \
				'specified by taxonomy keywords bounded by quotes: for AND use ' \
				'character \'&\', for OR use character \'|\', for NOT use character ' \
				'\'-\'. (E.g. "Actinobacteria|Proteobacteria-Pseudomonadales" which ' \
				'means Actinobacteria OR Proteobacteria without (NOT) ' \
				'Pseudomonadales). Please note that in order to limit OFFLINE BLASTP search ' \
				'by taxomony, user needs to provide file with taxids (see manual). [default: "Bacteria"]',
			type=str, required=False, default="Bacteria",
			metavar="\"BLASTP_TAXONOMY\"",
			)
		parser.add_argument('-k','--ktaxonomy',
			help='Specified genera/species (genera is represented by just one ' \
				'word; species is represented by just two words separated by space) ' \
				'will be represented by 3 peptides (per genera/species) in a final ' \
				'list of selected peptides. Bound specified taxonomy by quotes; for ' \
				'AND use character \'&\' (E.g. "Pseudomonas&Pseudomonas putida" will ' \
				'output 3 peptides representing genera Pseudomonas in file ' \
				'\'Selected_peptides.genera\' as well as 3 peptides representing ' \
				'species Pseudomonas putida in file \'Selected_peptides.species\'). ' \
				'[default: None]',
			type=str, required=False, default=None, metavar="\"KEY_TAXONOMY\"",
			)
		parser.add_argument('-s', '--specificity',
			help='Peptide specificity score threshold. [default: 90.0]',
			type=float, required=False, default=90.0, metavar="SPECIFICITY_THRESHOLD",
			)
		parser.add_argument('-a', '--avoid-spec-filt-list',
			help='Provide a list of peptides for that specificity scores are omitted. ' \
				'I.e. specificity score-based filtering step will not apply for such peptides. ' \
				'Bound specified peptides by quotes and separate individual entries by comma ' \
				'(E.g. "ENPPVLPK,SGLFTSEELPR").',
			type=str, required=False, default=None,
			)
		parser.add_argument('-cd', '--cdetectability',
			help='Minimum number of CONSeQence ML algorithms (out of 4) considering a peptide ' \
				'as detectable by MS to pass detectability filter. Please note that maximum ' \
				'value is 4. [default: 1]',
			type=int, required=False, default=1,
			)
		parser.add_argument('-pd', '--pdetectability',
			help='PeptideSieve detectability score threshold. [default: 0.6]',
			type=float, required=False, default=0.6, metavar="DETECTABILITY_THRESHOLD",
			)
		parser.add_argument('-f', '--chemmod-filter',
			help='Minimum number of possible chemical modifications of a peptide ' \
				'to be omitted from final selection lists. [default: 4]',
			type=int, required=False, default=4,
			)
		parser.add_argument('-n', '--peptide-calc-count',
			help='Number of peptides with the highest coverage in terms of microbial species ' \
				'to calculate specificity and detectability scores for. [default: 400]',
			type=int, required=False, default=400,
			)
		parser.add_argument('-N', '--peptide-out-count',
			help='Number of peptides to be output to \'Selected_peptides.*\' files. [default: 50]',
			type=int, required=False, default=50,
			)
		parser.add_argument('-t', '--threads',
			help='Number of threads (processes) to be used in parallel blasting. ' \
				'Works for online as well as offline blasting. [default: 4]',
			type=int, required=False, default=4,
			)
		parser.add_argument('-T', '--timeout',
			help='Maximum time (in seconds) to wait for online BLASTP results. ' \
				'[default: 1500]',
			type=int, required=False, default=1500,
			)
		parser.add_argument('-z', '--min-synname-length',
			help='Minimum number of words blastp synonymous name consist of to be ' \
				'accepted as synonymous name. [default: 1]',
			type=int, required=False, default=1,
			)
		parser.add_argument('-U', '--in-url',
			help='Specify UniProt URL in quotes. Applicable only if --custom-url ' \
				'argument is specified.',
			type=str, required='--custom-url' in sys.argv, default=None,
			)
		parser.add_argument('-i', '--in-aa-file',
			help='Specify path to dataset with protein sequences in FASTA AA format (e.g. ' \
				'infile.faa). Applicable only if --no-usearch argument is specified. Following  ' \
                'header formats are allowed: >PROTEIN_ID or >PROTEIN_ID~COV=<value> if ' \
                'coverage (expected copy number) of entries is also provided. ' \
                'Value can be either integer or float larger than 1.0.',
			type=str, required='--no-usearch' in sys.argv, default=None,
			)
		parser.add_argument('-e', '--in-extra-aa-file',
			help='Specify path to an extra input fasta AA file (e.g. infile_extra.faa). ' \
				'Applicable only if --extra-input argument is specified. Following  ' \
                'header formats are allowed: >PROTEIN_ID or >PROTEIN_ID~COV=<value> in ' \
                'case of coverage (expected copy number) of entries is also provided. ' \
                'Value can be either integer or float larger than 1.0.',
			type=str, required='--extra-input' in sys.argv, default=None,
			)
		parser.add_argument('-I', '--in-blastp-file',
			help='Specify path to \'blastp.xml\' file generated in previous program ' \
				'run. Applicable only if --no-bsearch argument is specified.',
			type=str, required='--no-bsearch' in sys.argv, default=None,
			)
		parser.add_argument('--stats-only',
			help='Print general info, protein statistics and figures, and exit.',
			required=False, action='store_true', default=False,
			)
		return parser.parse_args()

	def get_resultsFolder_name(self, mol):
		return '%s/results_%s_%s' % (pathlib.Path(sys.argv[0]).parent.absolute(),
			re.sub('[^a-zA-Z0-9]', '', mol).lower(),
			strftime('%d_%b_%Y_%H.%M.%S', localtime()))

	def parse_config_file(self, configPath):
		parser = ConfigParser()
		cfgExists = os.path.isfile(configPath)
		if cfgExists:
			parser.read(configPath)
		else:
			sys.stdout.write('Caution! Config file does not exist on a specified' \
			 	' path \'%s\'!\n' % configPath)
			sys.stdout.flush()
		return parser

	def check_for_data_in_config(self, configData, section, keyword):
		if configData.has_section(section):
			if configData.has_option(section, keyword):
				return True
			else:
				sys.stdout.write('Caution! Option %s is missing in section %s in a config file!\n' % (keyword, section))
				sys.stdout.flush()
				return False
		else:
			sys.stdout.write('Caution! Section %s is missing in a config file!\n' % section)
			sys.stdout.flush()
			return False

	def get_uniprot_url(self, molecule, taxonomy):
		url = 'https://www.uniprot.org/uniprot/?query=name%3A'
		isFirst = True
		for entry in re.split('(?<=\D)\+', re.sub('[^a-zA-Z0-9]+', '+', molecule)):
			if isFirst:
				if '+' in entry:
					url += '\"%s\"' % entry
				else:
					url += '%s' % entry
				isFirst = False
			else:
				if '+' in entry:
					url += '+name%%3A%%22%s%%22' % entry
				else:
					url += '+name%%3A%s' % entry
		for entry in re.split('([&|-])', taxonomy):
			if str(entry) == '&':
				url += '+AND'
			elif str(entry) == '|':
				url += '+OR'
			elif str(entry) == '-':
				url += '+NOT'
			else:
				url += '+taxonomy%%3A%s' % str(entry)
		url += '&columns=id,protein%20names,organism,lineage(PHYLUM),length,sequence'
		return url

	def retrieve_raw_uniprot_input(self, url):
		sys.stdout.write('Retrieving raw input from UniProt database...\n')
		sys.stdout.flush()
		try:
			rawData = urllib.request.urlopen(
				urllib.request.Request(url, urllib.parse.urlencode({'format':'tab'}).encode("utf-8"))).read().decode("utf-8")
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Cannot connect to UniProt server. Exiting... ' \
				'Try again later.\n')
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return rawData

	def extract_uniprot_proteins(self, rawData, mol, noPutatives):
		sys.stdout.write('Extracting %s entries...\n' % mol)
		sys.stdout.flush()
		linePattern = re.compile(r'(.+)\n')
		uniprotTabPattern = re.compile(r'(.+?)\t(.+?)\t(.+?)\t(.*?)\t(.+?)\t(.+)')
		for line in linePattern.findall(rawData):
			data = uniprotTabPattern.findall(line)
			if data:
				if data[0][0] != 'Entry':
					if noPutatives:
						if ("PUTATIVE" not in str(data[0][1]).upper() and
								"PROBABLE" not in str(data[0][1]).upper()):
							if (not data[0][3] or
									str(data[0][3].split(' ')[0]).upper() == 'CANDIDATE'):
								self.proteinEntries.append([
									data[0][0], 1.0, data[0][2], 'Unknown',
									data[0][4], data[0][5], '',
									])
							else:
								self.proteinEntries.append([
									data[0][0], 1.0, data[0][2], data[0][3],
									data[0][4], data[0][5], '',
									])
					else:
						if (not data[0][3] or
								str(data[0][3].split(' ')[0]).upper() == 'CANDIDATE'):
							self.proteinEntries.append([
								data[0][0], 1.0, data[0][2], 'Unknown',
								data[0][4], data[0][5], '',
								])
						else:
							self.proteinEntries.append([
								data[0][0], 1.0, data[0][2], data[0][3],
								data[0][4], data[0][5], '',
								])
			else:
				sys.stderr.write('Error: Wrong UniProt tab data format. Probably some' \
					' Uniprot data are missing. Exiting...\n')
				return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		if not self.proteinEntries:
			sys.stderr.write('Error: Raw UniProt data were obtained, however,' \
				' proteinEntries variable is empty. Exiting...\n')
			return False
		return True

	def get_uniprot_data(self, rawData, mol, noPutatives):
		if not rawData:
			sys.stderr.write('Error: URL request returned no data. Probable protein' \
				' name or taxonomy specification problem. Exiting...\n')
			return False
		if not self.extract_uniprot_proteins(rawData, mol, noPutatives):
			return False
		return True

	def print_extracted_uniprot_proteins(self):
		sys.stdout.write('Printing proteins extracted from the raw UniProt data...\n')
		sys.stdout.flush()
		try:
			with open('%s/UniProt_proteins.faa' % self.resultsFolder, 'w') as protFile:
				try:
					for protein in self.proteinEntries:
						protFile.write('>%s|%s|%s\n%s\n' % (protein[0], protein[2], protein[3], protein[5]))
				except:
					traceback.print_exc()
					sys.stderr.write('Error: Unexpected error occurred while writing proteins into file %s. Exiting...\n' % protFile)
					return False
		except OSError:
			sys.stderr.write('Error: Could not open %s file for writing. Exiting...\n' % protFile)
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_fasta_data(self, infile):
		sys.stdout.write('Extracting FASTA protein entries...\n')
		sys.stdout.flush()
		try:
			inFile = open(infile, 'r')
			try:
				headerPatternWithCov = re.compile(r'^>(.+?)~COV=(\d+(?:\.\d+)?)\n')
				headerPatternSimple = re.compile(r'^>(.+?)\n')
				header = []
				sequence = ''
				headerData = headerPatternWithCov.findall(inFile.readline())
				if headerData:
					if float(headerData[0][1]) < 1.0:
						sys.stderr.write('Error: %s entry has coverage lower than 1.0.' \
						' Exitting..' % headerData[0][0])
						return False
					header = [str(headerData[0][0]), float(headerData[0][1])]
				else:
					inFile.seek(0)
					headerData = headerPatternSimple.findall(inFile.readline())
					if headerData:
						header = [str(headerData[0]), 1.0]
					else:
						sys.stderr.write('Error: %s file has wrong format. Only fasta' \
    						' amino acid format (*.faa) is accepted. Following header' \
                            ' formats are allowed: >PROTEIN_ID or >PROTEIN_ID~COV=<value> in' \
                            ' case of coverage of entries is also provided. Value can be' \
                            ' either integer or float larger than 1.0. Exiting...\n' % infile)
						return False
				for line in inFile:
					headerData = headerPatternWithCov.findall(line)
					if headerData:
						if float(headerData[0][1]) < 1.0:
							sys.stderr.write('Error: %s entry has coverage lower than 1.0.' \
							' Exitting..' % line)
							return False
						headerNew = [str(headerData[0]), float(headerData[0][1])]
						if sequence:
							self.proteinEntries.append([
								header[0], header[1], 'Unknown', 'Unknown', len(sequence), sequence, '',
								])
							header = headerNew
							sequence = ''
						else:
							sys.stderr.write('Error: Missing AA sequence in a file %s.' \
								'Exiting...\n' % infile)
							return False
					else:
						headerData = headerPatternSimple.findall(line)
						if headerData:
							headerNew = [str(headerData[0]), 1.0]
							if sequence:
								self.proteinEntries.append([
									header[0], header[1], 'Unknown', 'Unknown', len(sequence), sequence, '',
									])
								header = headerNew
								sequence = ''
							else:
								sys.stderr.write('Error: Missing AA sequence in a file %s.' \
									'Exiting...\n' % infile)
								return False
						else:
							sequence += line.rstrip()
				self.proteinEntries.append([
					header[0], header[1], 'Unknown', 'Unknown', len(sequence), sequence, '',
					])
			except:
				traceback.print_exc()
				sys.stderr.write('Error: Unexpected error occurred while parsing file %s. Exiting...\n' % infile)
				inFile.close()
				return False
		except OSError:
			sys.stderr.write('Error: Could not open %s file. Exiting...\n' % infile)
			return False
		finally:
			inFile.close()
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_extra_data(self, infile):
		sys.stdout.write('Extracting CDS entries from extra FASTA file...\n')
		extraProteins = []# sequence, coverage
		try:
			inFile = open(infile, 'r')#
			try:
				headerPatternWithCov = re.compile(r'^>(.+?)~COV=(\d+(?:\.\d+)?)\n')
				headerPatternSimple = re.compile(r'^>(.+?)\n')
				protCov = 0.0
				sequence = ''
				headerData = headerPatternWithCov.findall(inFile.readline())
				if headerData:
					if float(headerData[0][1]) < 1.0:
						sys.stderr.write('Error: %s entry has coverage lower than 1.0.' \
						' Exitting..' % headerData[0][0])
						return False
					protCov = float(headerData[0][1])
				else:
					inFile.seek(0)
					headerData = headerPatternSimple.findall(inFile.readline())
					if headerData:
						protCov = 1.0
					else:
						sys.stderr.write('Error: %s file has wrong format. Only fasta' \
    						' amino acid format (*.faa) is accepted. Following header' \
                            ' formats are allowed: >PROTEIN_ID or >PROTEIN_ID~COV=<value> in' \
                            ' case of coverage of entries is also provided. Value can be' \
                            ' either integer or float larger than 1.0. Exiting...\n' % infile)
						return False
				for line in inFile:
					headerData = headerPatternWithCov.findall(line)
					if headerData:
						if float(headerData[0][1]) < 1.0:
							sys.stderr.write('Error: %s entry has coverage lower than 1.0.' \
							' Exitting..' % line)
							return False
						protCovNew = float(headerData[0][1])
						if sequence:
							extraProteins.append([sequence, protCov])
							protCov = protCovNew
							sequence = ''
						else:
							sys.stderr.write('Error: Missing AA sequence in a file %s.' \
								'Exiting...\n' % infile)
							return False
					else:
						headerData = headerPatternSimple.findall(line)
						if headerData:
							protCovNew = 1.0
							if sequence:
								extraProteins.append([sequence, protCov])
								protCov = protCovNew
								sequence = ''
							else:
								sys.stderr.write('Error: Missing AA sequence in a file %s.' \
									'Exiting...\n' % infile)
								return False
						else:
							sequence += line.rstrip()
				extraProteins.append([sequence, protCov])
			except:
				traceback.print_exc()
				sys.stderr.write('Error: Unexpected error occurred while parsing file %s. Exiting...\n' % infile)
				inFile.close()
				return False
		except OSError:#
			sys.stderr.write('Error: Could not open %s file. Exiting...\n' % infile)
			return False
		finally:
			inFile.close()
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return extraProteins

	def get_input_data(self, args):
		os.makedirs(self.resultsFolder)
		if args.no_usearch:
			if not self.get_fasta_data(args.in_aa_file):
				return False
		elif args.custom_url:
			customUrl = args.in_url + '&columns=id,protein%20names,organism,lineage(PHYLUM),length,sequence'
			if not self.get_uniprot_data(
					self.retrieve_raw_uniprot_input(customUrl),
					args.molecule, args.no_putatives):
				return False
		else:
			if not self.get_uniprot_data(
					self.retrieve_raw_uniprot_input(
					self.get_uniprot_url(args.molecule, args.utaxonomy)),
					args.molecule, args.no_putatives):
				return False
		return True

	def draw_len_histogram(self, title, minLength, maxLength, removedProteinLengths):
		try:
			if maxLength - minLength > 1500:
				divider = 3
			else:
				divider = 2
			if removedProteinLengths:
				if ((removedProteinLengths[-1] - removedProteinLengths[0]) / divider) > 0:
					remNumOfBins = int((removedProteinLengths[-1] - removedProteinLengths[0])
									/ divider)
				else:
					remNumOfBins = 1
			keptNumOfBins = int((int(self.proteinEntries[-1][4]) - int(self.proteinEntries[0][4]))
							/ divider)
			if keptNumOfBins == 0:
				keptNumOfBins = 1
			fig, ax = matplotlib.pyplot.subplots()
			ax.hist(
				[int(item[4]) for item in self.proteinEntries],
				color='blue',
				bins=keptNumOfBins,
				label='kept',
				)
			if removedProteinLengths:
				ax.hist(removedProteinLengths, color='red', bins=remNumOfBins, label='discarded')
			ax.set(xlabel='length [AA]', ylabel='count')
			matplotlib.pyplot.title(title, fontsize=10, fontweight='bold')
			matplotlib.pyplot.legend(loc='upper right')
			matplotlib.pyplot.savefig('%s/length_distribution.uniprot_entries.png' % self.resultsFolder, dpi=400)
			return True
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Unexpected error occurred while drawing histogram. Exiting...\n')
			return False

	def draw_pie_chart(self, title, labels, fractions, figName, fontSize, labelDist):
		matplotlib.pyplot.figure(figsize=(6,6))
		matplotlib.rcParams.update({'font.size': fontSize})
		matplotlib.pyplot.title(title, fontsize=10, fontweight='bold')
		_, _, autotexts = matplotlib.pyplot.pie(
			fractions,
			labels=labels,
			autopct='%1.1f%%',
			labeldistance=labelDist,
			colors=[
				'blue', 'green', 'red', 'goldenrod',
				'teal', 'magenta', 'olive', 'black',
				'orangered', 'purple', 'grey', 'sienna',
				])
		for autotext in autotexts:
			autotext.set_color('white')
		matplotlib.pyplot.savefig('%s/' % self.resultsFolder + figName, dpi=400)

	def extract_lower_taxonomy(self, proteinClass, inEntries):
		unknownGeneraSynonyms = ['BLOOD', 'UNCULTURED', 'CANDIDATE', 'BACTERIUM', 'BACTERIA']# invalid uniprot genera entries
		unknownSpeciesSynonyms = ['sp.', 'genomosp.', 'bacterium', 'bacteria']# invalid uniprot species entries
		try:
			if str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']') == 'Candidatus':
				genusAlreadyContained = False
				for genus in inEntries[1:]:
					if (genus[0] == '%s %s' % (
							str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']'),
							str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']'))):
						genusAlreadyContained = True
						if (str(proteinClass[2].split(' ')[2]).lstrip('[').rstrip(']') not in
								unknownSpeciesSynonyms):
							if (str(proteinClass[2].split(' ')[2]).lstrip('[').rstrip(']') not in
									genus[1]):
								genus[1][str(proteinClass[2].split(' ')[2]).lstrip('[').rstrip(']')] = 1
							else:
								genus[1][str(proteinClass[2].split(' ')[2]).lstrip('[').rstrip(']')] += 1
						else:
							if 'unknown' not in genus[1]:
								genus[1]['unknown'] = 1
							else:
								genus[1]['unknown'] += 1
						break
				if not genusAlreadyContained:
					if (str(proteinClass[2].split(' ')[2]).lstrip('[').rstrip(']') not in
							unknownSpeciesSynonyms):
						inEntries.append([
							'%s %s' % (
								str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']'),
								str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']')),
							{str(proteinClass[2].split(' ')[2]).lstrip('[').rstrip(']'): 1},
							])
					else:
						inEntries.append([
							'%s %s' % (
								str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']'),
								str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']')),
							{'unknown': 1},
							])
			elif (str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']').upper() in
				unknownGeneraSynonyms):
				genusAlreadyContained = False
				for genus in inEntries[1:]:
					if genus[0] == 'Unknown':
						genus[1]['unknown'] += 1
						genusAlreadyContained = True
						break
				if not genusAlreadyContained:
					inEntries.append(['Unknown', {'unknown': 1}])
			else:
				genusAlreadyContained = False
				for genus in inEntries[1:]:
					if genus[0] == str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']'):
						genusAlreadyContained = True
						if (str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']') not in
								unknownSpeciesSynonyms):
							if (str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']') not in
									genus[1]):
								genus[1][str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']')] = 1
							else:
								genus[1][str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']')] += 1
						else:
							if 'unknown' not in genus[1]:
								genus[1]['unknown'] = 1
							else:
								genus[1]['unknown'] += 1
						break
				if not genusAlreadyContained:
					if (str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']') not in
							unknownSpeciesSynonyms):
						inEntries.append([
							str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']'),
							{str(proteinClass[2].split(' ')[1]).lstrip('[').rstrip(']'): 1},
							])
					else:
						inEntries.append([
							str(proteinClass[2].split(' ')[0]).lstrip('[').rstrip(']'),
							{'unknown': 1},
							])
		except:
			genusAlreadyContained = False
			for genus in inEntries[1:]:
				if genus[0] == 'Unknown':
					genus[1]['unknown'] += 1
					genusAlreadyContained = True
					break
			if not genusAlreadyContained:
				inEntries.append(['Unknown', {'unknown': 1}])

	def extract_taxonomy_stats(self):
		sys.stdout.write('Extracting taxonomy statistics...\n')
		sys.stdout.flush()
		try:
			for phylum in Counter([item[3] for item in self.proteinEntries]).most_common():
				self.phylumEntries.append([phylum[0]])
				for protein in self.proteinEntries:
					if protein[3] == phylum[0]:
						self.extract_lower_taxonomy(protein, self.phylumEntries[-1])
			for phylum in self.phylumEntries:
				for genus in phylum[1:]:
					genus[1] = OrderedDict(sorted(list(genus[1].items()), key=lambda x: x[1], reverse=True))
				phylum[1:] = sorted(
					phylum[1:],
					key=lambda x_dic: sum(x_dic[1].values()),
					reverse=True)
			sys.stdout.write('Done.\n\n')
			sys.stdout.flush()
			return True
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Unexpected error occurred while extracting taxonomy stats. Exiting...\n')
			return False

	def print_general_stats(
			self, args, allCount, allMedian,
			allAvg, allMin, allMax, allStd,
			removedProteinLengths):
		sys.stdout.write('Printing general statistics and drawing figures...\n')
		sys.stdout.flush()
		try:
			with open('%s/initial_info.txt' % self.resultsFolder, 'w') as statsFile:
				statsFile.write('Arguments: --no-usearch %s\n' \
				 				'           --extra-input %s\n' \
								'           --custom-url %s\n' \
								'           --no-bsearch %s\n' \
					            '           --offline-blastp %s\n' \
								'           --run-peptidesieve %s\n' \
								'           -m "%s"\n' \
								'           -C %s\n' \
								'           -l %.3f\n' \
					            '           --no-putatives %s\n' \
								'           -u "%s"\n' \
								'           -c %d\n' \
								'           -x %d\n' \
								'           -y %d\n' \
								'           -b "%s"\n' \
								'           -k "%s"\n' \
								'           -s %.3f\n' \
								'           -a %s\n' \
								'           -cd %d\n' \
								'           -pd %.3f\n' \
					            '           -f %d\n' \
								'           -n %d\n' \
								'           -N %d\n' \
								'           -t %d\n' \
								'           -T %d\n' \
								'           -z %d\n' \
								'           -U %s\n' \
								'           -i %s\n' \
								'           -e %s\n' \
								'           -I %s\n' \
								'           --stats-only %s\n\n' % (
						args.no_usearch, args.extra_input, args.custom_url, args.no_bsearch,
						args.offline_blastp, args.run_peptidesieve, args.molecule, args.config,
						args.length, args.no_putatives, args.utaxonomy, args.cleavages,
						args.mins, args.maxs, args.btaxonomy, args.ktaxonomy,
						args.specificity, args.avoid_spec_filt_list, args.cdetectability, args.pdetectability,
						args.chemmod_filter, args.peptide_calc_count, args.peptide_out_count, args.threads,
						args.timeout, args.min_synname_length, args.in_url, args.in_aa_file,
						args.in_extra_aa_file, args.in_blastp_file,args.stats_only)
						)
				statsFile.write('### INITIAL UNIPROT/CUSTOM DATA ###\n\n')
				statsFile.write('Total number of proteins:\t%d\nMedian length [AA]:\t%.1f\n' \
					'Average length [AA]:\t%.1f\nMin length [AA]:\t%d\nMax length [AA]:\t%d\n' \
					'Standard deviation [AA]:\t%.2f\n' % (
					allCount, allMedian, allAvg, allMin, allMax, allStd))
				statsFile.write('\n### REFINED UNIPROT/CUSTOM DATA ###\n\n')
				statsFile.write('Total number of proteins:\t%d\nMedian length [AA]:\t%.1f\n' \
					'Average length [AA]:\t%.1f\nMin length [AA]:\t%d\nMax length [AA]:\t%d\n' \
					'Standard deviation [AA]:\t%.2f\n' % (
						len(self.proteinEntries),
						numpy.median([int(item[4]) for item in self.proteinEntries]),
						numpy.mean([int(item[4]) for item in self.proteinEntries]),
						int(self.proteinEntries[0][4]),
						int(self.proteinEntries[-1][4]),
						numpy.std([int(item[4]) for item in self.proteinEntries])))
				statsFile.write('\n### TAXONOMY OF REFINED DATA ###\n\n')
				statsFile.write('Protein counts:\n\n')
				phylaLabels = []
				phylaFractions = []
				fractOther = 0
				phylaCount = 0
				for phylum in self.phylumEntries:
					phylaCount += 1
					spCount = 0
					for genus in phylum[1:]:
						spCount += sum(genus[1].values())
					statsFile.write('%s phylum: %s\n\nGenera:\n\n' % (phylum[0].upper(), spCount))
					if spCount*100 / len(self.proteinEntries) > 3.0:
						phylaLabels.append(phylum[0])
						phylaFractions.append(spCount*100 / len(self.proteinEntries))
					else:
						fractOther += spCount
					generaLabels = []
					generaFractions = []
					fractOtherG = 0
					for genus in phylum[1:]:
						statsFile.write('%s: %d (species:' % (genus[0], sum(genus[1].values())))
						for key, value in genus[1].items():
							statsFile.write(' %s: %d;' % (key, value))
						statsFile.write(')\n')
						if phylaCount < 3:
							if sum(genus[1].values())*100 / spCount > 1.5:
								generaLabels.append(genus[0])
								generaFractions.append(sum(genus[1].values())*100 / spCount)
							else:
								fractOtherG += sum(genus[1].values())
					if fractOtherG > 0:
						generaLabels.append('Other')
						generaFractions.append(fractOtherG*100 / spCount)
					if generaFractions:
						self.draw_pie_chart(
							phylum[0],
							generaLabels,
							generaFractions,
							'genera_pie_chart_%d.uniprot_entries.png' % phylaCount,
							6,
							1.03)
					statsFile.write('\n\n')
				if fractOther > 0:
					phylaLabels.append('Other')
					phylaFractions.append(fractOther*100/len(self.proteinEntries))
				self.draw_pie_chart(
					args.molecule,
					phylaLabels,
					phylaFractions,
					'phyla_pie_chart.uniprot_entries.png',
					10,
					1.1)
				if not self.draw_len_histogram(args.molecule, allMin, allMax, removedProteinLengths):
					return False
				sys.stdout.write('Done.\n\n')
				sys.stdout.flush()
				return True
		except OSError:
			sys.stderr.write('Error: Could not open %s/initial_info.txt file for writing.' % self.resultsFolder)
			return False
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Unexpected error occurred while printing general stats. Exiting...\n')
			return False

	def preprocess_proteins(self, args):
		sys.stdout.write('Filtering proteins according to their length...\n')
		sys.stdout.flush()
		discardedProteins = []
		self.proteinEntries.sort(key=lambda x: int(x[4]))
		count = len(self.proteinEntries)
		median = numpy.median([int(item[4]) for item in self.proteinEntries])
		mean = numpy.mean([int(item[4]) for item in self.proteinEntries])
		std = numpy.std([int(item[4]) for item in self.proteinEntries])
		minL = int(self.proteinEntries[0][4])
		maxL = int(self.proteinEntries[-1][4])
		for entry in self.proteinEntries:
			if (int(entry[4]) > args.length * median + median or
					int(entry[4]) < median - args.length * median):
				discardedProteins.append(entry)
		for entry in discardedProteins:
			self.proteinEntries.remove(entry)
		if not self.proteinEntries:
			sys.stderr.write('Error: All protein entries were removed during initial filtering steps. ' \
				'Try to modify argument values. Exiting...\n')
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		if not self.extract_taxonomy_stats():
			return False
		if not self.print_general_stats(
				args, count, median, mean, minL, maxL, std,
				[int(item[4]) for item in discardedProteins]):
			return False
		return True

	def digest_protein(self, seq, cleav, minPep, maxPep):
		rawPeptides = []
		resPeptides = []
		try:
			aaSeq = AASequence.fromString(seq)
			digestor = ProteaseDigestion()
			digestor.setMissedCleavages(cleav)
			digestor.digest(aaSeq, rawPeptides)
		except:
			sys.stdout.write('Caution! There are problems with the digestion of protein %s.' \
				'Omitting this protein...\n' % seq)
			sys.stdout.flush()
			return ''
		for entry in rawPeptides:
			if len(entry.toString()) >= minPep and len(entry.toString()) <= maxPep:
				resPeptides.append(entry.toString())
		return resPeptides

	def run_protein_digestion(self, cleav, minPep, maxPep):
		sys.stdout.write('Protein digestion in progress...\n')
		sys.stdout.flush()
		for entry in self.proteinEntries:
			entry[6] = self.digest_protein(entry[5], cleav, minPep, maxPep)
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()

	def run_extra_protein_digestion(self, cleav, minPep, maxPep):
		sys.stdout.write('Extra protein dataset digestion in progress...\n')
		sys.stdout.flush()
		extraPeptides = []# each entry consist of aa sequence, frequency, different species count, peptideSieve detectability score, CONSeQuence detectability, specificity score, possible chemical modifications count
		for entry in self.extraProteinEntries:
			for peptide in self.digest_protein(entry[0], cleav, minPep, maxPep):
				isInEPlist = False
				for subEPlist in extraPeptides:
					if subEPlist[0] == peptide and len(subEPlist[0]) == len(peptide):
						subEPlist[1] += int(round(entry[1]))
						isInEPlist = True
						break
				if not isInEPlist:
					extraPeptides.append([peptide.decode(), int(round(entry[1])), 1, None, None, None, None])
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return extraPeptides

	def calculate_peptide_frequencies(self):
		return Counter(list(chain(*[item[6] for item in self.proteinEntries]))).most_common()

	def extract_peptide_info(self, noUsearch):
		sys.stdout.write('Extracting peptide info...\n')
		sys.stdout.flush()
		try:
			for peptide in self.calculate_peptide_frequencies():
				if 'Windows' in platform.system():
					self.peptideEntries.append([peptide[0], peptide[1], None, None, None, None, None]) # sequence, frequency, different species count, peptideSieve detectability score, CONSeQuence detectability, specificity score, possible chemical modifications count
				else:
					self.peptideEntries.append([peptide[0].decode(), peptide[1], None, None, None, None, None]) # sequence, frequency, different species count, peptideSieve detectability score, CONSeQuence detectability, specificity score, possible chemical modifications count
				for protein in self.proteinEntries:
					if protein[6].count(peptide[0]) > 0:
						phylumAlreadyContained = False
						if protein[1] > 1.0:
							self.peptideEntries[-1][1] += int(round(protein[6].count(peptide[0])*(protein[1]-1.0),0))
						for phylum in self.peptideEntries[-1][7:]:
								if str(protein[3]) == phylum[0]:
									self.extract_lower_taxonomy(protein, phylum)
									phylumAlreadyContained = True
						if not phylumAlreadyContained:
							self.peptideEntries[-1].append([str(protein[3])])
							self.extract_lower_taxonomy(protein, self.peptideEntries[-1][-1])
				uniqueSpeciesCount = 0
				for phylum in self.peptideEntries[-1][7:]:
					for genus in phylum[1:]:
						uniqueSpeciesCount += len(genus[1])
				self.peptideEntries[-1][2] = uniqueSpeciesCount
			if noUsearch:
				self.peptideEntries.sort(key=lambda x: x[1], reverse=True)
			self.peptideEntries.sort(key=lambda x: x[2], reverse=True)
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Unexpected error occured. Exiting...\n')
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_peptideSieve_detectability_values(self, peptideSieveFile):
		sys.stdout.write('Extracting peptideSieve detectability values...\n')
		sys.stdout.flush()
		detectValPattern = re.compile(r'.+?\t.+?\t(.+?)\t(.+?)\n')
		try:
			with open(peptideSieveFile, 'r') as inFile:
				for line in inFile:
					data = detectValPattern.findall(line)
					if data:
						for peptide in self.peptideEntries:
							if str(data[0][0]) == peptide[0]:
								peptide[3] = float(data[0][1])
								break
					else:
						sys.stderr.write('Error: Could not parse %s file.' % peptideSieveFile)
						return False
		except OSError:
			sys.stderr.write('Error: Could not open %s file.' % peptideSieveFile)
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_relevant_peptides(self):
		sys.stdout.write('Extracting relevant peptides...\n')
		sys.stdout.flush()
		relevantPeptides = []
		for extraPeptide in self.extraPeptideEntries:
			for peptide in self.peptideEntries:
				if extraPeptide[0] == peptide[0] and len(extraPeptide[0]) == len(peptide[0]):
					relevantPeptides.append(extraPeptide)
		relevantPeptides.sort(key=lambda x: x[1], reverse=True)
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return relevantPeptides

	def print_get_most_frequent_relevant_peptides(self, countLimit, splitIntoMultiFiles):
		sys.stdout.write('Printing most relevant peptides to \*txt and \*faa files...\n')
		sys.stdout.flush()
		try:
			with open('%s/peptide_candidates_for_calculations.txt' % self.resultsFolder, 'w') as txtFile:
				peptideCount = 1
				fileCount = 1
				peptCountPerFile = 1
				try:
					faFile = open('%s/peptide_candidates_for_calculations%.2d.faa' % (
						self.resultsFolder, fileCount), 'w')
					for relPeptide in self.relevantPeptides:
						if 'X' in str(relPeptide[0]):
							sys.stdout.write('Caution: Omitting peptide %s since it contains amino acid X.\n' % str(relPeptide[0]))
							sys.stdout.flush()
							continue
						if splitIntoMultiFiles:
							if peptCountPerFile % 51 == 0:
								fileCount += 1
								peptCountPerFile = 1
								faFile.close()
								try:
									faFile = open('%s/peptide_candidates_for_calculations%.2d.faa' % (
										self.resultsFolder, fileCount), 'w')
								except OSError:
									sys.stderr.write('Error: Could not open %s/peptide_candidates_for_calculations%.2d.faa ' \
										'file for writing. Exiting...\n' % (self.resultsFolder, fileCount))
									return False
						txtFile.write('peptide_%d_freq_count_%d\t%s\n' % (
							peptideCount, int(relPeptide[1]), str(relPeptide[0])))
						faFile.write('>peptide_%d_freq_count_%d\n%s\n' % (
							peptideCount, int(relPeptide[1]), str(relPeptide[0])))
						self.queryPeptides.append(['>peptide_%d_freq_count_%d' % (
							peptideCount, int(relPeptide[1])), str(relPeptide[0])])
						peptideCount += 1
						peptCountPerFile += 1
						if peptideCount > countLimit:
							faFile.close()
							break
					faFile.close()
					sys.stdout.write('Done.\n\n')
					sys.stdout.flush()
					return True
				except OSError:
					sys.stderr.write('Error: Could not open %s/peptide_candidates_for_calculations%.2d.faa file for writing. ' \
						'Exiting...\n' % (self.resultsFolder, fileCount))
					return False
		except OSError:
			sys.stderr.write('Error: Could not open %s/peptide_candidates_for_calculations.txt file for writing. Exiting...\n'
				% self.resultsFolder)
			return False

	def print_get_most_frequent_peptides(self, countLimit, splitIntoMultiFiles):
		sys.stdout.write('Printing most frequent peptides to *txt and *faa files...\n')
		sys.stdout.flush()
		try:
			with open('%s/peptide_candidates_for_calculations.txt' % self.resultsFolder, 'w') as txtFile:
				peptideCount = 1
				fileCount = 1
				peptCountPerFile = 1
				try:
					faFile = open('%s/peptide_candidates_for_calculations%.2d.faa' % (
						self.resultsFolder, fileCount), 'w')
					for peptide in self.peptideEntries:
						if 'X' in str(peptide[0]):
							sys.stdout.write('Caution: Omitting peptide %s since it contains amino acid X.\n' % str(peptide[0]))
							sys.stdout.flush()
							continue
						if splitIntoMultiFiles:
							if peptCountPerFile % 51 == 0:
								fileCount += 1
								peptCountPerFile = 1
								faFile.close()
								try:
									faFile = open('%s/peptide_candidates_for_calculations%.2d.faa' % (
										self.resultsFolder, fileCount), 'w')
								except OSError:
									sys.stderr.write('Error: Could not open %s/peptide_candidates_for_calculations%.2d.faa file for ' \
										'writing. Exiting...\n' % (self.resultsFolder, fileCount))
									return False
						txtFile.write('peptide_%d_diff_sp_count_%d\t%s\n' % (
							peptideCount, int(peptide[2]), str(peptide[0])))
						faFile.write('>peptide_%d_diff_sp_count_%d\n%s\n' % (
							peptideCount, int(peptide[2]), str(peptide[0])))
						self.queryPeptides.append(['>peptide_%d_diff_sp_count_%d' % (
							peptideCount, int(peptide[2])), str(peptide[0])])
						peptideCount += 1
						peptCountPerFile += 1
						if peptideCount > countLimit:
							faFile.close()
							break
					faFile.close()
					sys.stdout.write('Done.\n\n')
					sys.stdout.flush()
					return True
				except OSError:
					sys.stderr.write('Error: Could not open %s/peptide_candidates_for_calculations%.2d.faa file for writing. ' \
						'Exiting...\n' % (self.resultsFolder, fileCount))
					return False
		except OSError:
			sys.stderr.write('Error: Could not open %s/peptide_candidates_for_calculations.txt file for writing. Exiting...\n'
				% self.resultsFolder)
			return False

	def run_consequence(self, isExtraInput, driverPath):
		sys.stdout.write('Running CONSeQuence...\n')
		sys.stdout.flush()
		try:
			options = webdriver.ChromeOptions()
			options.add_argument("headless")
			sys.stdout.write('Launching Google Chrome browser in the background.\n')
			sys.stdout.flush()
			driver = webdriver.Chrome(executable_path=driverPath, options=options)
		except:
			if 'Darwin' in platform.system():
				try:
					sys.stdout.write('Failed to launch Google Chrome browser. Is your chromedriver compatible with Chrome?\n')
					sys.stdout.write('Darwin-based OS detected. Launching Safari browser. Caution: Cannot be ran in the background.\n')
					sys.stdout.flush()
					driver = webdriver.Safari()
				except:
					traceback.print_exc()
					sys.stderr.write('Failed to launch Safari browser.\n')
					sys.stderr.write('Error: Could not open any of web browsers. Please check the manual for troubleshooting.')
					return False
			else:
				traceback.print_exc()
				sys.stderr.write('Error: Failed to launch Google Chrome browser. Is your chromedriver compatible with Chrome? Please ' \
					'check the manual for troubleshooting.')
				return False
		try:
			buffer = ""
			for qPept in self.queryPeptides:
				buffer += "%s\n%s\n" % (qPept[0], qPept[1])
			driver.get('http://king.smith.man.ac.uk/CONSeQuence/')
			assert "CONSeQuence" in driver.title
			pageSource = driver.page_source
			elem = driver.find_element_by_name("sequence")
			elem.clear()
			elem.send_keys(buffer)
			elem.submit()
			while True:
				time.sleep(5)
				pageSource = driver.page_source
				status = re.search(r'div class=\"main\"', pageSource)
				if status:
					break
			resPattern = re.compile(r'.+?PEPTIDE.+?</td><td>(.+?)</td><td>0</td><td>(\d)</td>')
			resPeptides = resPattern.findall(driver.page_source)
			if not resPeptides:
				sys.stderr.write('Error: CONSeQuence did not return the results. It provided following output instead (in HTML format):\n')
				sys.stderr.write(driver.page_source)
				return False
			if isExtraInput:
				peptideEntries = self.relevantPeptides
			else:
				peptideEntries = self.peptideEntries
			for qPept in self.queryPeptides:
				isInResPeptides = False
				for rPept in resPeptides:
					if rPept[0] == qPept[1]:
						for peptide in peptideEntries:
							if peptide[0] == qPept[1]:
								peptide[4] = int(rPept[1])
								isInResPeptides = True
								break
				if not isInResPeptides:
					for peptide in peptideEntries:
						if peptide[0] == qPept[1]:
							peptide[4] = 0
							break
			sys.stdout.write('Done.\n\n')
			sys.stdout.flush()
			return True
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Unexpected error occurred while accessing CONSeQuence server programmatically.\n')
			return False

	def run_peptideSieve(self, inFile, executable, propertiesPath):
		sys.stdout.write('Running PeptideSieve...\n')
		sys.stdout.flush()
		try:
			subprocess.call([
				executable,
				'-e', '.peptideSieveOUT',
				'-f', 'TXT',
				'-d', 'MUDPIT_ESI',
				'-p', '0.0',
				'-P', propertiesPath,
				'-O', '%s' % self.resultsFolder,
				inFile])
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Could not run PeptideSieve.\n')
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def detect_possible_chem_modifications(self, isExtraInput, countLimit):
		sys.stdout.write('Detecting peptides prone to chemical modifications...\n')
		sys.stdout.flush()
		try:
			with open('%s/possible_chemical_modifications.txt' % self.resultsFolder, 'w') as chmFile:
				chmFile.write('Peptide\t#\tPossible chemical modifications\n')
				peptCount = 0
				peptideEntries = []
				if isExtraInput:
					for peptide in self.relevantPeptides:
						peptideEntries.append(peptide)
						peptCount += 1
						if peptCount == countLimit:
							break
				else:
					for peptide in self.peptideEntries:
						peptideEntries.append(peptide)
						peptCount += 1
						if peptCount == countLimit:
							break
				for peptide in peptideEntries:
					chmCount = 0
					chmDescript = ''
					if peptide[0][0] == 'Q' or peptide[0][0] == 'E':
						chmCount += 1
						if chmDescript == '':
							chmDescript += 'pyroE'
						else:
							chmDescript += '; pyroE'
					if peptide[0][0] == 'C':
						if chmDescript == '':
							chmDescript += 'pyro-cmC'
						else:
							chmDescript += '; pyro-cmC'
					problematicSubstrs = ['W', 'M', 'C', 'DP', 'DG', 'NG', 'QG']
					for substr in problematicSubstrs:
						if peptide[0].count(substr) != 0:
							chmCount += peptide[0].count(substr)
							if chmDescript == '':
								chmDescript += '%s(%d)' % (substr, peptide[0].count(substr))
							else:
								chmDescript += '; %s(%d)' %(substr, peptide[0].count(substr))
					chmFile.write('%s\t%d\t%s\n' % (peptide[0], chmCount, chmDescript))
					peptide[6] = chmCount
		except OSError:
			sys.stderr.write('Error: Could not open %s/possible_chemical_modifications.txt file for writing. Exiting...\n'
				% self.resultsFolder)
			return False
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Unexpected error occurred. Exiting...\n')
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_entrezQuery(self, taxonomy):
		entrezQuery = ""
		for entry in re.split('([&|-])', taxonomy):
			if str(entry) == '&':
				entrezQuery += " AND "
			elif str(entry) == '|':
				entrezQuery += " OR "
			elif str(entry) == '-':
				entrezQuery += " NOT "
			else:
				entrezQuery += "%s[ORGN]" % str(entry)
		return entrezQuery

	def run_blastp(self, executable, databasePath, inFile, threads):
		sys.stdout.write('Offline blastp search in progress...\n')
		sys.stdout.flush()
		try:
			subprocess.call([
				executable,
				'-db', databasePath,
				'-evalue', '1000000.0',
				'-max_target_seqs', '100',
				'-word_size', '4',
				'-gapopen', '9',
				'-gapextend', '1',
				'-num_threads', str(threads),
				'-matrix', 'PAM30',
				'-threshold', '11',
				'-outfmt', '5',
				'-taxidlist', '/Users/matejmedvecky/Documents/hcmr/pipeline/2to3/out/marinomonas.ids',
				'-query', inFile,
				'-out', '%s/blastp_results.xml' % self.resultsFolder])
		except:
			traceback.print_exc()
			sys.stderr.write('Error: Could not run offline blastp program.\n')
			return False
		sys.stdout.write('Blastp search is done.\n')
		sys.stdout.flush()
		return True

	def run_qblast(self, program, database, inFile, entrezQuery):
		reconnectCount = 0
		while True:
			try:
				resHandle = NCBIWWW.qblast(
					program,
					database,
					open(inFile).read(),
					expect=1000000.0,
					hitlist_size=100,
					word_size=2,
					gapcosts="9 1",
					matrix_name='PAM30',
					threshold=11,
					entrez_query=entrezQuery)
				try:
					f = open('%s/blastp_results.xml' % self.resultsFolder, 'a')
					f.write(resHandle.read())
					f.close()
					resHandle.close()
					return True
				except OSError:
					sys.stderr.write('Error: Could not open %s/blastp_results.xml file for appending. Exiting...\n'
						% self.resultsFolder)
					return False
				except:
					traceback.print_exc()
					sys.stderr.write('Error: Unexpected error occurred while writing results into file blastp_results.xml.' \
						'Exiting...\n')
					return False
			except:
				sys.stdout.write('Cannot connect to NCBI server, trying to reconnect...\n')
				sys.stdout.flush()
				time.sleep(15)
				reconnectCount += 1
				if reconnectCount == 4:
					sys.stderr.write('Error: There are problems with the connection to NCBI' \
						' server. Blastp analysis of peptides could not be accomplished.\n')
					return False

	def perform_qblast_search(self, program, database, inFiles, entrezQuery, threads, timeout):
		sys.stdout.write('Online blastp search in progress...\n')
		sys.stdout.flush()
		resendCount = 0
		while resendCount < 3:
			if resendCount != 0:
				sys.stdout.write('Waiting for the response from NCBI server for too long,' \
					' resending the queries...\n')
				sys.stdout.flush()
			qblastSearchPool = multiprocess.Pool(threads)
			multi_res = [qblastSearchPool.apply_async(
				self.run_qblast,
				args = (program, database, inFile, entrezQuery)) for inFile in inFiles]
			for res in multi_res:
				toBeResended = False
				try:
					res.get(timeout)
				except multiprocess.TimeoutError:
					qblastSearchPool.terminate()
					qblastSearchPool.join()
					resendCount += 1
					toBeResended = True
					if resendCount == 3:
						sys.stderr.write('Warning: There are problems obtaining results from' \
							' NCBI sever. At least some of them are missing.\n')
						sys.stderr.flush()
						toBeResended = False
					break
			if toBeResended:
				try:
					os.remove('%s/blastp_results.xml' % self.resultsFolder)
				except OSError:
					pass
				continue
			break
		sys.stdout.write('Blastp search is done.\n')
		sys.stdout.flush()
		return True

	def refine_blast_file(self, inBlastFile):#get rid of lines without XML-like pattern which raise ValueError during file parsing
		sys.stdout.write('Refining blast XML file...\n')
		sys.stdout.flush()
		try:
			with open(inBlastFile, 'r') as inFile:
				try:
					with open('%s/blastp_results_refined.xml' % self.resultsFolder, 'w') as outFile:
						xmlPattern = re.compile(r'.*?<.+?>\n')
						for line in inFile:
							xmlData = xmlPattern.findall(line)
							if xmlData:
								outFile.write(line)
				except OSError:
					sys.stderr.write('Error: Could not open %s/blastp_results_refined.xml file for writing. Exiting...\n' % self.resultsFolder)
					return False
		except OSError:
			sys.stderr.write('Error: Could not open %s file for reading. Exiting...\n' % inBlastFile)
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_peptide_specificity(self, mol, inBlastFile, minSynNameLength, isExtraInput):
		sys.stdout.write('Computing peptide specificity...\n')
		sys.stdout.flush()
		try:
			with open(inBlastFile, 'r') as inFile:
				try:
					with open('%s/peptide_blastp_specificity.txt' % self.resultsFolder, 'w') as specFile:
						synNames = [re.sub('[^a-zA-Z0-9]+', ' ', mol).upper(),]
						currentPept = []
						if isExtraInput:
							peptideEntries = self.relevantPeptides
						else:
							peptideEntries = self.peptideEntries
						peptMissingInBlastResFile = []
						alnMissingInBlastResFile = []
						for record in NCBIXML.parse(inFile):
							if not record.alignments:
								alnMissingInBlastResFile.append(record.query)
								continue
							isInBlastResults = False
							for qPept in self.queryPeptides:
								if qPept[0].lstrip('>') == record.query:
									currentPept = qPept
									isInBlastResults = True
							if not isInBlastResults:
								peptMissingInBlastResFile.append(record.query)
								continue
							for align in record.alignments:
								for hsp in align.hsps:
									if hsp.match == currentPept[1] and record.query_length == len(currentPept[1]):
										if (re.search(
												re.sub('[^a-zA-Z0-9]', '', mol).upper(),
												re.sub('[^a-zA-Z0-9]', '', align.title).upper())):
											for name in re.split('>.+?\s|;|(?<=\D), ', align.title):
												name = re.sub(
													'\[.+?\]|gi\|.+?\s|ref\|.+?\s|gb\|.+?\s|emb\|.+?\s|tpg\|.+?\s' \
														'|dbj\|.+?\s|sp\|.+?\s|MULTISPECIES:|RecName: Full=|AltName: Full=',
													'',
													name)
												if (re.sub('[^a-zA-Z0-9]+', ' ', name).upper().strip() not in
														synNames and not re.search(
														'HYPOTHETICAL PROTEIN|UNCHARACTERIZED PROTEIN|' \
															'UNCHARACTERIZED CONSERVED PROTEIN|PARTIAL',
														re.sub('[^a-zA-Z0-9]+', ' ', name).upper().strip())):
													wordCount = 0
													isSynonymous = False
													for word in re.split(
															'(?<=\D) ',
															re.sub('[^a-zA-Z0-9]+', ' ', name).upper().strip()):
														wordCount += 1
														if re.search(r'\b%s\b' % re.escape(word), re.sub(
																'[^a-zA-Z0-9]+', ' ', mol).upper()):
															isSynonymous = True
													if wordCount >= minSynNameLength and isSynonymous:
														synNames.append(re.sub(
															'[^a-zA-Z0-9]+', ' ', name).upper().strip())
									break##
						if alnMissingInBlastResFile:
							sys.stdout.write('Warning: Following peptides are lacking alignments in' \
								' blastp_results.xml file, therefore they are omitted from' \
								' specificity score calculations. Try to re-run the analysis' \
								' in order to get the alignments.\n%s\n\n' % alnMissingInBlastResFile)
							sys.stdout.flush()
						if peptMissingInBlastResFile:
							sys.stderr.write('Caution: Following list of query peptide headers from provided' \
								' blastp_results.xml file do not occur within the list of peptide entries' \
								' digested and selected from input proteins. Did you specify correct' \
								' blastp_results.xml file?\n%s\n\n' % peptMissingInBlastResFile)
							sys.stdout.flush()
						specFile.write('List of detected synonymous names to \'%s\':\n' % mol)
						for synName in synNames:
							specFile.write('%s\n' % synName.lower())
						specFile.write('\n# # # # # # # # # # # # # # # # # # # # # # # # # # #' \
							' # # # # # # # # # # # # # # # # # # # # # #\n\n')
						inFile.seek(0)
						currentPept = []
						for record in NCBIXML.parse(inFile):
							fullMatchCount = 0
							goodProteinCount = 0
							hypoProteinCount = 0
							specScore = [100.0, 0.0]
							diffNames = []
							if not record.alignments:
								continue
							isInBlastResults = False
							for qPept in self.queryPeptides:
								if qPept[0].lstrip('>') == record.query:
									currentPept = qPept
									isInBlastResults = True
							if not isInBlastResults:
								continue
							for align in record.alignments:
								matchedSynName = False
								for hsp in align.hsps:
									if hsp.match == currentPept[1] and record.query_length == len(currentPept[1]):
										fullMatchCount += 1
										for name in synNames:
											if (re.search(name, re.sub(
													'[^a-zA-Z0-9]+', ' ', align.title).upper())):
												goodProteinCount += 1
												matchedSynName = True
												break
										if not matchedSynName:
											if re.search(
													'HYPOTHETICAL PROTEIN|UNCHARACTERIZED PROTEIN|' \
														'UNCHARACTERIZED CONSERVED PROTEIN',
													re.sub('[^a-zA-Z0-9]', ' ', align.title).upper()):
												hypoProteinCount += 1
											else:
												diffNames.append(align.title)
									break
							specFile.write('Peptide: %s\nNumber of full query sequence matches: %d/100\n' % (
								currentPept[1], fullMatchCount))
							if fullMatchCount != 0:
								specScore[1] = 100.0 / fullMatchCount
							if diffNames:
								bestMatches = []
								for diffName in diffNames:
									matches = []
									bestMatches.append(['', '', 0, 0])# best diffName, synName, shared words count, number of synName's words
									for name in re.split('>.+?\s|;|(?<=\D), ', diffName):
										name = re.sub(
											'\[.+?\]|gi\|.+?\s|ref\|.+?\s|gb\|.+?\s|emb\|.+?\s|tpg\|.+?\s|dbj\|.+?\s|sp\|.+?\s' \
											'|MULTISPECIES:|RecName: Full=|AltName: Full=',
											'',
											name)
										for synName in synNames:
											matches.append([
												re.sub('[^a-zA-Z0-9]+', ' ', name).upper().strip(),
												synName,
												0,
												0,
												])
											for entry in re.split('(?<=\D) ', synName):
												matches[-1][3] += 1
											for word in re.split(
													'(?<=\D) ',
													re.sub('[^a-zA-Z0-9]+', ' ', name).upper().strip()):
												if re.search(r'\b%s\b' % re.escape(word), synName):
													matches[-1][2] += 1
										for match in matches:
											if match[2] == 0 and bestMatches[-1][2] == 0:
												if bestMatches[-1][0] == '':
													bestMatches[-1][0] = match[0]
												else:
													nameAlreadyContained = False
													for entry in re.split(' \| ', bestMatches[-1][0]):
														if entry == match[0]:
															nameAlreadyContained = True
													if not nameAlreadyContained:
														bestMatches[-1][0] += ' | ' + match[0]
											elif match[2] != 0 and bestMatches[-1][2] == 0:
												bestMatches[-1][0] = match[0]
												bestMatches[-1][1] = match[1]
												bestMatches[-1][2] = match[2]
												bestMatches[-1][3] = match[3]
											elif (((float(match[2]) / match[3])
												> (float(bestMatches[-1][2]) / bestMatches[-1][3])) or
												((float(match[2]) / match[3])
												== (float(bestMatches[-1][2]) / bestMatches[-1][3]) and
												match[2] > bestMatches[-1][2])):
												bestMatches[-1][0] = match[0]
												bestMatches[-1][1] = match[1]
												bestMatches[-1][2] = match[2]
												bestMatches[-1][3] = match[3]
								bestMatches.sort(key = lambda x: int(x[2]), reverse = True)
								diffProteinCount = 0
								partDiffProteinCount = 0
								for bMatch in bestMatches:
									if bMatch[2] > 0:
										if bMatch[2] == bMatch[3]:
											goodProteinCount += 1
										else:
											partDiffProteinCount += 1
									else:
										diffProteinCount += 1
								specFile.write('Number of \'%s\' or synonymous protein names in' \
									' full query matches: %d\nNumber of hypothetical proteins in' \
									' full query matches: %d\n' % (
									mol, goodProteinCount, hypoProteinCount))
								specFile.write('Number of partially different protein names in' \
									' full query matches: %d\nNumber of different protein names in' \
									' full query matches: %d\n' % (
									partDiffProteinCount, diffProteinCount))
								specScore[0] -= diffProteinCount * specScore[1]
								if bestMatches[0][2] > 0:
									firstPartDiffProt = True
									for bMatch in bestMatches:
										if bMatch[2] > 0:
											if bMatch[2] != bMatch[3]:
												if firstPartDiffProt == True:
													firstPartDiffProt = False
													specFile.write('List of partially different' \
														' protein names:\n')
												specFile.write('%s --> matches %d/%d words in %s\n' % (
													bMatch[0], bMatch[2], bMatch[3], bMatch[1]))
												specScore[0] -= ((float(bMatch[2])
																/ bMatch[3])
																* specScore[1])
								if bestMatches[-1][2] == 0:
									specFile.write('List of different protein names:\n')
									for bMatch in bestMatches:
										if bMatch[2] == 0:
											specFile.write('%s\n' % (bMatch[0]))
							else:
								specFile.write('Number of \'%s\' or synonymous protein names in' \
									' full query matches: %d\nNumber of hypothetical proteins in' \
									' full query matches: %d\n' % (
									mol, goodProteinCount, hypoProteinCount))
								specFile.write('Number of partially different protein names in' \
									' full query matches: 0\nNumber of different protein names in' \
									' full query matches: 0\n')
							specFile.write('Score out of 100: %.1f\n\n' % specScore[0])
							isInPeptideEntries = False
							for entry in peptideEntries:
								if entry[0] == currentPept[1] and len(currentPept[1]) == len(entry[0]):
									entry[5] = specScore[0]
									isInPeptideEntries = True
									break
							if not isInPeptideEntries:
								sys.stderr.write('Caution: Query peptide %s from provided' \
								 	' blastp_results.xml file does not occur within the list' \
									' of peptide entries generated from input proteins. Did' \
									' you specify correct blastp_results.xml file?\n' % hsp.query)
								sys.stdout.flush()
				except OSError:
					sys.stderr.write('Error: Could not open %s/peptide_blastp_specificity.txt file for writing. Exiting...\n'
						% self.resultsFolder)
					return False
				except:
					traceback.print_exc()
					sys.stderr.write('Error: Unspecified error while parsing file %s and calculating '
						'specificity values. Exiting...\n' % inFile)
					return False
		except OSError:
			sys.stderr.write('Error: Could not open %s file for reading. Exiting...\n' % inBlastFile)
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def print_peptides_complete_info(self, isExtraInput):
		sys.stdout.write('Printing peptide info...\n')
		sys.stdout.flush()
		try:
			with open('%s/Peptides_complete_info.txt' % self.resultsFolder, 'w') as peptFile:
				peptFile.write('sequence\tfrequency\tdiff sp count\tpeptideSieve detectability score\t' \
					'CONSeQuence detectability\tspecificity score\tpossible chemical modifications count\t' \
					'representing taxonomy\n')
				if isExtraInput:
					peptideEntries = self.relevantPeptides
					sys.stdout.write('Caution: Metagenomic input. Omitting peptides derived from Uniprot ' \
						'that do not match metagenomic data.\n')
					sys.stdout.flush()
				else:
					peptideEntries = self.peptideEntries
				for peptide in peptideEntries:
					peptFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (
						peptide[0], peptide[1], peptide[2], peptide[3], peptide[4], peptide[5], peptide[6]))
					if not peptide[7:]:
						peptFile.write('N/A\n')
					else:
						for phylum in peptide[7:]:
							peptFile.write('%s: ' % (phylum[0]))
							for genus in phylum[1:]:
								peptFile.write('%s' % (genus[0]))
								isFirst = True
								for key in genus[1]:
									if isFirst:
										peptFile.write(' %s' % key)
										isFirst = False
									else:
										peptFile.write(', %s' % key)
								peptFile.write('; ')
							peptFile.write('\t')
						peptFile.write('\n')
		except OSError:
			sys.stderr.write('Error: Could not open %s/Peptides_complete_info.txt file for writing. Exiting...\n'
				% self.resultsFolder)
			return False
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

	def get_print_highest_coverage_peptides(
			self, getGeneraCoverage, isExtraInput, isNoUsearch, outCount, keyTaxonomy,
			specificityThres, avoidSpecFiltString, pDetectabilityThres, cDetectabilityThres,
			ChModCount):
		keyTax = []
		if getGeneraCoverage:
			sys.stdout.write('Getting highest genera coverage peptides...\n')
			sys.stdout.flush()
			if keyTaxonomy:
				for entry in keyTaxonomy.split('&'):
					if ' ' not in entry:
						keyTax.append([entry, 0])
		else:
			sys.stdout.write('Getting highest species coverage peptides...\n')
			sys.stdout.flush()
			if keyTaxonomy:
				for entry in keyTaxonomy.split('&'):
					if ' ' in entry:
						keyTax.append([entry, 0])
		listOfNames = []#genus/species, is covered
		selectedPeptides = []
		bestNextCandidate = ['', []]#peptide, added genera/species entries
		avoidSpecFiltList = []
		if avoidSpecFiltString:
			sys.stdout.write('Caution: Specificity score-based filtering step will not apply to following list ' \
				'of peptides specified as a parameter: %s.\n' % avoidSpecFiltString)
			sys.stdout.flush()
			avoidSpecFiltList = avoidSpecFiltString.split(",")
		for phylum in self.phylumEntries:
			for genus in phylum[1:]:
				if genus[0] != 'Unknown':
					if getGeneraCoverage and genus[0] not in [item[0] for item in listOfNames]:#get rid of potential taxonomy duplications (e.g. due to assignment of genera / species to wrong phyla in some UniProt DB entries)
						listOfNames.append([genus[0], False])
					else:
						for key, value in genus[1].items():
							if '%s %s' % (genus[0], key) not in [item[0] for item in listOfNames]:#get rid of potential taxonomy duplications (e.g. due to assignment of genera / species to wrong phyla in some UniProt DB entries)
								listOfNames.append(['%s %s' % (genus[0], key), False])
		for entry in keyTax:
			if entry[0] not in [item[0] for item in listOfNames]:
				sys.stdout.write('Caution: There is no peptide representing \'%s\'' \
					' specified as key taxonomy member.\n' % entry[0])
				sys.stdout.flush()
				keyTax.remove(entry)
		if isExtraInput:
			peptideEntries = self.relevantPeptides
		else:
			peptideEntries = self.peptideEntries
		unsuitablePeptCount = 0
		for peptide in peptideEntries:
			if 'X' in peptide[0]:
				unsuitablePeptCount += 1
				continue
			if peptide[3] != None:
				if peptide[3] < pDetectabilityThres:
					unsuitablePeptCount += 1
					continue
			if peptide[4] != None:
				if peptide[4] < cDetectabilityThres:
					unsuitablePeptCount += 1
					continue
			if peptide[5] != None:
				isInAvoidList = False
				if avoidSpecFiltList:
					for avoidSpecFiltPeptide in avoidSpecFiltList:
						if (str(avoidSpecFiltPeptide).strip() == peptide[0] and
							len(str(avoidSpecFiltPeptide).strip()) == len(peptide[0])):
								isInAvoidList = True
								break
				if not isInAvoidList:
					if peptide[5] < specificityThres:
						unsuitablePeptCount += 1
						continue
			if peptide[6] != None:
				if peptide[6] >= ChModCount:
					unsuitablePeptCount += 1
					continue
		if (len(peptideEntries) - unsuitablePeptCount) < outCount:
			sys.stdout.write('Caution: There is not enough peptides passing filtering' \
				' conditions. --peptide-out-count value was changed to %d.\n' % (len(peptideEntries) - unsuitablePeptCount))
			sys.stdout.flush()
			outCount = (len(peptideEntries) - unsuitablePeptCount)
		while len(selectedPeptides) < outCount:
			if keyTax:
				noMoreKeyPeptides = True
			for peptide in peptideEntries:
				coversKeyTaxonomy = False
				if (peptide[0] in [item[0] for item in selectedPeptides] or
					peptide[2] <= len(bestNextCandidate[1])):
						continue
				if 'X' in peptide[0]:
					continue
				if peptide[3] != None:
					if peptide[3] < pDetectabilityThres:
						continue
				if peptide[4] != None:
					if peptide[4] < cDetectabilityThres:
						continue
				if peptide[5] != None:
					isInAvoidList = False
					if avoidSpecFiltList:
						for avoidSpecFiltPeptide in avoidSpecFiltList:
							if (str(avoidSpecFiltPeptide).strip() == peptide[0] and
								len(str(avoidSpecFiltPeptide).strip()) == len(peptide[0])):
									isInAvoidList = True
									break
					if not isInAvoidList:
						if peptide[5] < specificityThres:
							continue
				if peptide[6] != None:
					if peptide[6] >= ChModCount:
						continue
				addedNames = []
				if keyTax:
					for phylumPept in peptide[7:]:
						for genusPept in phylumPept[1:]:
							if getGeneraCoverage:
								if genusPept[0] in [genusKey[0] for genusKey in keyTax]:
									noMoreKeyPeptides = False
									coversKeyTaxonomy = True
								for genusEntry in listOfNames:
									if (genusPept[0] == genusEntry[0] and not genusEntry[1] and
											genusPept[0] not in addedNames):
										addedNames.append(genusEntry[0])
										break
							else:
								for key, value in genusPept[1].items():
									if ('%s %s' % (genusPept[0], key) in
												[speciesKey[0] for speciesKey in keyTax]):
										noMoreKeyPeptides = False
										coversKeyTaxonomy = True
									for speciesEntry in listOfNames:
										if ('%s %s' % (genusPept[0], key) == speciesEntry[0] and not
												speciesEntry[1] and '%s %s' % (genusPept[0], key) not in
												addedNames):
											addedNames.append(speciesEntry[0])
											break
					if not coversKeyTaxonomy:
						continue
					else:
						if not bestNextCandidate[0]:
							bestNextCandidate[0] = peptide[0]
							bestNextCandidate[1] = addedNames
						elif len(bestNextCandidate[1]) < len(addedNames):
							bestNextCandidate[0] = peptide[0]
							bestNextCandidate[1] = addedNames
				else:
					for phylumPept in peptide[7:]:
						for genusPept in phylumPept[1:]:
							if getGeneraCoverage:
								for genusEntry in listOfNames:
									if (genusPept[0] == genusEntry[0] and not genusEntry[1] and
											genusPept[0] not in addedNames):
										addedNames.append(genusEntry[0])
										break
							else:
								for key, value in genusPept[1].items():
									for speciesEntry in listOfNames:
										if ('%s %s' % (genusPept[0], key) == speciesEntry[0] and not
												speciesEntry[1] and '%s %s' % (genusPept[0], key) not in
												addedNames):
											addedNames.append(speciesEntry[0])
											break
					if not bestNextCandidate[0]:
						bestNextCandidate[0] = peptide[0]
						bestNextCandidate[1] = addedNames
					elif len(bestNextCandidate[1]) < len(addedNames):
						bestNextCandidate[0] = peptide[0]
						bestNextCandidate[1] = addedNames
			if keyTax:
				if noMoreKeyPeptides:
					for keyTaxon in keyTax:
						sys.stderr.write('Caution: %d suitable peptide/s only was/were' \
							' found for %s taxon.\n' % (keyTaxon[1], keyTaxon[0]))
						sys.stderr.flush()
						keyTax.remove(keyTaxon)
				else:
					for peptide in peptideEntries:
						if peptide[0] == bestNextCandidate[0]:
							for phylumPept in peptide[7:]:
								for genusPept in phylumPept[1:]:
									if getGeneraCoverage:
										for genusKey in keyTax:
											if genusKey[0] == genusPept[0]:
												genusKey[1] += 1
												if genusKey[1] == 3:
													keyTax.remove(genusKey)
									else:
										for key, value in genusPept[1].items():
											for speciesKey in keyTax:
												if speciesKey[0] == '%s %s' % (genusPept[0], key):
													speciesKey[1] += 1
													if speciesKey[1] == 3:
														keyTax.remove(speciesKey)
			if bestNextCandidate[0] != '':
				selectedPeptides.append(bestNextCandidate)
				for addedName in bestNextCandidate[1]:
					for nameEntry in listOfNames:
						if addedName == nameEntry[0]:
							nameEntry[1] = True
							break
			bestNextCandidate = ['', []]
		if isExtraInput or isNoUsearch:
			outFileName = 'Selected_peptides.no_taxonomy'
		elif getGeneraCoverage:
			outFileName = 'Selected_peptides.genera'
		else:
			outFileName = 'Selected_peptides.species'
		try:
			outFile = open('%s/%s' % (self.resultsFolder, outFileName), 'w')
		except OSError:
			sys.stderr.write('Error: Could not open %s/%s file for writing. Exiting...\n' % (self.resultsFolder, outFileName))
			return False
		if not isExtraInput and not isNoUsearch:
			if getGeneraCoverage:
				if listOfNames:
					outFile.write('Genera coverage: %.1f %%\n\nSelected peptides/added genera:\n' % (
						(float([genus[1] for genus in listOfNames].count(True)) / len(listOfNames))*100))
				else:
					outFile.write('Genera coverage: N/A\n\nSelected peptides/added genera:\n')
			else:
				if listOfNames:
					outFile.write('Species coverage: %.1f %%\n\nSelected peptides/added species:\n' % (
						(float([species[1] for species in listOfNames].count(True)) / len(listOfNames))*100))
				else:
					outFile.write('Species coverage: N/A\n\nSelected peptides/added species:\n')
		else:
			outFile.write('Selected peptides:\n')
		selPeptWithNoSpecScores = []
		selPeptWithNoDetScores = []
		for peptide in selectedPeptides:
			outFile.write('%s\t' % peptide[0])
			for entry in peptideEntries:
				if entry[0] == peptide[0] and len(entry[0]) == len(peptide[0]):
					if entry[5] == None:
						selPeptWithNoSpecScores.append(peptide[0])
					if entry[4] == None:
						selPeptWithNoDetScores.append(peptide[0])
						sys.stderr.flush()
					break
			isFirst = True
			for species in peptide[1]:
				if isFirst:
					outFile.write(species)
					isFirst = False
				else:
					outFile.write(', ' + species)
			outFile.write('\n')
		outFile.close()
		if selPeptWithNoSpecScores:
			sys.stderr.write('Caution: Following list of selected peptides is missing specificity scores' \
				' (they are probably out of --peptide-calc-count range). Therefore they' \
				' could be nonspecific!\n%s\n' % selPeptWithNoSpecScores)
			sys.stderr.flush()
		if selPeptWithNoDetScores:
			sys.stderr.write('Caution: Following list of selected peptides is missing CONSeQuence' \
				' detectability values (they are probably out of --peptide-calc-count' \
				' range). Therefore they could be not amenable to LC-MS/MS analysis!\n%s\n' % selPeptWithNoDetScores)
			sys.stderr.flush()
		sys.stdout.write('Done.\n\n')
		sys.stdout.flush()
		return True

def main():
	mp = MetaproteomicPipline()
	mp.arguments = mp.parse_user_input()
	mp.resultsFolder = mp.get_resultsFolder_name(mp.arguments.molecule)
	mp.configData = mp.parse_config_file(mp.arguments.config)
	if not mp.get_input_data(mp.arguments):
		return False
	if not mp.arguments.no_usearch:
		if not mp.print_extracted_uniprot_proteins():
			return False
	if not mp.preprocess_proteins(mp.arguments):
		return False
	if not mp.arguments.stats_only:
		mp.run_protein_digestion(mp.arguments.cleavages, mp.arguments.mins, mp.arguments.maxs)
		if not mp.extract_peptide_info(mp.arguments.no_usearch):
			return False
		if mp.arguments.extra_input:
			mp.extraProteinEntries = mp.get_extra_data(mp.arguments.in_extra_aa_file)
			if not mp.extraProteinEntries:
				return False
			mp.extraPeptideEntries = mp.run_extra_protein_digestion(mp.arguments.cleavages, mp.arguments.mins, mp.arguments.maxs)
			mp.relevantPeptides = mp.get_relevant_peptides()
			if not mp.relevantPeptides:
				sys.stderr.write('Error: None of peptides from metagenomic data match peptides derived from Uniprot database! Exitting.\n')
				return False
			if mp.arguments.offline_blastp:
				if not mp.print_get_most_frequent_relevant_peptides(mp.arguments.peptide_calc_count, False):
					return False
			else:
				if not mp.print_get_most_frequent_relevant_peptides(mp.arguments.peptide_calc_count, True):
					return False
		else:
			if mp.arguments.offline_blastp:
				if not mp.print_get_most_frequent_peptides(mp.arguments.peptide_calc_count, False):
					return False
			else:
				if not mp.print_get_most_frequent_peptides(mp.arguments.peptide_calc_count, True):
					return False
		if not mp.detect_possible_chem_modifications(mp.arguments.extra_input, mp.arguments.peptide_calc_count):
			return False
		if mp.check_for_data_in_config(mp.configData, 'ChromeDriver', 'path'):
			mp.chromedriverPath = mp.configData.get('ChromeDriver', 'path')
		else:
			mp.chromedriverPath = None
		if not mp.run_consequence(mp.arguments.extra_input, mp.chromedriverPath):
			return False
		if mp.arguments.run_peptidesieve:
			if mp.check_for_data_in_config(mp.configData, 'PeptideSieve', 'executable'):
				mp.ps_executable = mp.configData.get('PeptideSieve', 'executable')
				if mp.check_for_data_in_config(mp.configData, 'PeptideSieve', 'propertiesFilePath'):
					mp.propertiesFilePath = mp.configData.get('PeptideSieve', 'propertiesFilePath')
					if not mp.run_peptideSieve('%s/peptide_candidates_for_calculations.txt' % mp.resultsFolder, mp.ps_executable, mp.propertiesFilePath):
						return False
					if not mp.get_peptideSieve_detectability_values('%s/peptide_candidates_for_calculations.peptideSieveOUT' % mp.resultsFolder):
						return False
				else:
					sys.stderr.write('Therefore, PeptideSieve cannot be run. Skipping this step.\n')
					sys.stdout.flush()
			else:
				sys.stderr.write('Therefore, PeptideSieve cannot be run. Skipping this step.\n')
				sys.stdout.flush()
		if not mp.arguments.no_bsearch:
			if mp.arguments.offline_blastp:
				if mp.check_for_data_in_config(mp.configData, 'blastp', 'executable'):
					mp.blast_executable = mp.configData.get('blastp', 'executable')
					if mp.check_for_data_in_config(mp.configData, 'blastp', 'databasePath'):
						mp.blastDatabasePath = mp.configData.get('blastp', 'databasePath')
						if not mp.run_blastp(mp.blast_executable, mp.blastDatabasePath, '%s/peptide_candidates_for_calculations01.faa' % mp.resultsFolder, mp.arguments.threads):
							return False
					else:
						sys.stderr.write('Error: Could not run offline blastp program.\n')
						sys.stdout.flush()
						return False
				else:
					sys.stderr.write('Error: Could not run offline blastp program.\n')
					sys.stdout.flush()
					return False
			else:
				if not mp.perform_qblast_search('blastp', 'nr', [inFile for inFile in glob.glob('%s/*.faa' % mp.resultsFolder)], mp.get_entrezQuery(mp.arguments.btaxonomy), mp.arguments.threads, mp.arguments.timeout):
					return False
		if mp.arguments.no_bsearch and mp.arguments.in_blastp_file:
			if not mp.refine_blast_file(mp.arguments.in_blastp_file):
				return False
			if not mp.get_peptide_specificity(mp.arguments.molecule, '%s/blastp_results_refined.xml' % mp.resultsFolder, mp.arguments.min_synname_length, mp.arguments.extra_input):
				return False
		elif not mp.arguments.no_bsearch:
			if not mp.refine_blast_file('%s/blastp_results.xml' % mp.resultsFolder):
				return False
			if not mp.get_peptide_specificity(mp.arguments.molecule, '%s/blastp_results_refined.xml' % mp.resultsFolder, mp.arguments.min_synname_length, mp.arguments.extra_input):
				return False
		if not mp.print_peptides_complete_info(mp.arguments.extra_input):
			return False
		if not mp.get_print_highest_coverage_peptides(True, mp.arguments.extra_input, mp.arguments.no_usearch, mp.arguments.peptide_out_count, mp.arguments.ktaxonomy, mp.arguments.specificity, mp.arguments.avoid_spec_filt_list, mp.arguments.pdetectability, mp.arguments.cdetectability, mp.arguments.chemmod_filter):
			return False
		if not mp.arguments.extra_input or mp.arguments.no_usearch:
			if not mp.get_print_highest_coverage_peptides(False, mp.arguments.extra_input, mp.arguments.no_usearch, mp.arguments.peptide_out_count, mp.arguments.ktaxonomy, mp.arguments.specificity, mp.arguments.avoid_spec_filt_list, mp.arguments.pdetectability, mp.arguments.cdetectability, mp.arguments.chemmod_filter):
				return False
	sys.stdout.write('\nIf you use PepMANDIS in your research, please CITE our paper:\n\'Medvecky M, Mandalakis M. PepMANDIS: A peptide selection tool for designing function-based targeted proteomic assays in mixed microbial communities\'.\n\n')

if __name__ == '__main__':
	main()
