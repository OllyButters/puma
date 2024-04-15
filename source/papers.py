#!/usr/bin/env python3

################################################################################
# The publications metadata augmenter!
# This is the starting point of a pipeline that tries to augment a list of
# publications with metadata from places like PubMed and DOI.org to build a
# reporting tool and some pretty web pages.
# Go read the docs: https://github.com/OllyButters/puma/wiki
################################################################################

# core packages
import datetime
import json
import os.path
from os import listdir
import os
import logging
import time
import sys

# Internal packages
from setup import setup
from clean import clean
import add.geocode
import web_pages.build_htmlv2
import bibliography.bibtex
import get.simple_collate
from networks import author_network
from analyse import analyse
from analyse import coverage_report
from config import config

__author__ = "Olly Butters, Hugh Garner, Tom Burton, Becca Wilson"
__date__ = 12/4/2024
__version__ = '2.2'

# Time Log
start_time = time.time()
print('Start Time: ' + str(datetime.datetime.now().strftime("%H:%M")))

# Lets figure out some paths that everything is relative to
# global root_dir
path_to_papers_py = os.path.abspath(sys.argv[0])
root_dir = os.path.dirname(os.path.dirname(path_to_papers_py))
print('Root directory = ' + root_dir)

# Get all the config - these will be a global vars available like config.varname
config.build_config_variables(root_dir)

# Delete any unneeded data hanging around in the cache
setup.tidy_existing_file_tree()
setup.clean_old_zotero_cache_file()
setup.clean_old_pubmed_cache_file()
setup.clean_old_scopus_cache_file()

# Build the file tree relative to the root_dir
setup.build_file_tree()

# Set up the logging. Level is set in config and can be DEBUG, INFO, WARNING, ERROR, CRITICAL.
log_file = root_dir + '/logs/'+config.project_details['short_name']+'.log'
logging.basicConfig(level=config.logging_loglevel, filename=log_file,filemode="w",force=True)

print('Log file: ' + log_file)
print('Run something like: tail -f ' + log_file)

# Output some info to the log file to help with debugging
logging.info('Running version: %s', __version__)
logging.info('Started at: %s', str(datetime.datetime.now().strftime("%H:%M")))

###########################################################
# Get the metadata from external sources. This will store the raw metadata
# in a cache directory and the raw/merged metadata in the merged directory.
get.simple_collate.collate()


# Now read in all the cached merged metadata.
# 'papers' will be the giant LIST that has all the papers in it, each as a dictionary.
papers = []

# Get list of files in merged directory
merged_files_list = listdir(config.cache_dir + '/processed/merged/')
merged_files_list.sort()
# merged_files_list = merged_files_list[0:10]
print(str(len(merged_files_list))+' merged papers to load.')

# Open each one and add to papers object
for this_merged_file in merged_files_list:
    with open(config.cache_dir + '/processed/merged/' + this_merged_file, encoding='utf-8') as fo:
        # Will be a dictionary
        this_paper = json.load(fo)
        this_paper['filename'] = this_merged_file
        papers.append(this_paper)

print(str(len(papers)) + ' papers to process.')


###########################################################
# Clean the data - e.g. tidy dates and institute names
clean.clean(papers)

# should probably move clean_institution into main clean directly
clean.clean_institution(papers)

###########################################################
# Add some extra data in - i.e. geocodes and citations
add.geocode.geocode(papers)

# Write papers to summary file
file_name = root_dir + '/data/' + config.project_details['short_name'] + '/summary_added_to'
fo = open(file_name, 'w', encoding='utf-8')
fo.write(json.dumps(papers, indent=4))
fo.close()

# Write a copy of each paper to a separate file
for this_paper in papers:
    this_file_name = config.cache_dir + '/processed/cleaned/' + this_paper['IDs']['zotero'] + '.cleaned.json'
    fo = open(this_file_name, 'w', encoding='utf-8')
    fo.write(json.dumps(this_paper, indent=4))
    fo.close()

# Generate the coverage report.
coverage_report.coverage_report(papers)

# Output a BibTeX file with all the papers in it.
bibliography.bibtex.bibtex(papers)

###########################################################
# Do some actual analysis on the data. This will result in
# some CSV type files that can be analysed.
analyse.journals(papers)

# Figure out the word frequencies
analyse.word_frequencies(papers, 'title')
papers_with_keywords = analyse.word_frequencies(papers, 'keywords')
papers_with_abstract_text = analyse.word_frequencies(papers, 'abstract')

network = analyse.authors(papers)
analyse.first_authors(papers)
analyse.inst(papers)
# analyse.mesh(papers)
analyse.output_csv(papers)
analyse.dates(papers)

###########################################################
# Make some web pages
web_pages.build_htmlv2.build_css_colour_scheme()
age_weighted_citations, age_weighted_citations_data = web_pages.build_htmlv2.build_home(papers)
web_pages.build_htmlv2.build_papers(papers)
web_pages.build_htmlv2.build_mesh(papers)
web_pages.build_htmlv2.build_country_map(papers)
web_pages.build_htmlv2.build_metrics(papers, age_weighted_citations, age_weighted_citations_data, config.metrics_study_start_year, config.metrics_study_current_year)
web_pages.build_htmlv2.build_abstract_word_cloud(papers_with_abstract_text)
web_pages.build_htmlv2.build_keyword_word_cloud(papers_with_keywords)
web_pages.build_htmlv2.build_institute_map(papers)
web_pages.build_htmlv2.build_help()
web_pages.build_htmlv2.build_search(papers)


if config.web_page_show_zotero_tags:
    web_pages.build_htmlv2.build_zotero_tags(papers)

# Generate and dump the html for author network.
# Currently in development, not ready for general use.
if config.network_create_networks:
    web_pages.build_htmlv2.build_author_network(papers, network)
    author_network.build_network()

# Time Log
end_time = time.time()
elapsed_time = int(end_time - start_time)
print('End Time: ' + str(datetime.datetime.now().strftime("%H:%M")))
print('Elapsed Time (H:mm:ss) - ' + str(datetime.timedelta(seconds=elapsed_time)))

logging.info('End Time: %s', str(datetime.datetime.now().strftime("%H:%M")))
logging.info('Elapsed time (H:mm:ss) : %s', str(datetime.timedelta(seconds=elapsed_time)))
