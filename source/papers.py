#! /usr/bin/env python

import json
import os.path
from os import listdir
import logging
import time
# import pprint

# import get.get
import clean.clean
import add.geocode
import add.citations
import analyse.analysis
import html.htmlerrorlog.errorlog
import html.build_htmlv2
import bibliography.bibtex

##########################################################
# Get all the paper metadata from pubmed and do stuff with it
##########################################################
__author__ = "Olly Butters, Hugh Garner, Tom Burton"
__date__ = 27/7/16
__version__ = '0.2.6'


# Options - these should get moved out into a config file
# Stick a flag into see if we want to update the citations
update_citations = True
scopus_api_key = '8024d746590aade6be6856a22a734783'
scopus_citation_max_life = 30  # days
# pp = pprint.PrettyPrinter(indent=4)

# Time Log
start_time = time.time()
print "Start Time: " + str(start_time)


# Error log for displaying data input problems to user
error_log = html.htmlerrorlog.errorlog.ErrorLog()
# error_log.logError("Test Error")
# error_log.logWarning("Test Warning")


###########################################################
# Make sure the directory structure is set up first.
# Everything in the cache is grabbed from elsewhere, or built on the fly,
# so it should all be considerd ready to be deleted at any point!
if (os.path.exists('../cache') is False):
    os.mkdir('../cache')

if (os.path.exists('../data') is False):
    os.mkdir('../data')

# Log directory
if (os.path.exists('../logs') is False):
    os.mkdir('../logs')

# Output html
if not os.path.exists('../html'):
    os.mkdir('../html')

if not os.path.exists('../html/mesh'):
    os.mkdir('../html/mesh')

if not os.path.exists('../html/css'):
    os.mkdir('../html/css')

if not os.path.exists('../html/papers'):
    os.mkdir('../html/papers')

if not os.path.exists('../html/all_keywords'):
    os.mkdir('../html/all_keywords')

if not os.path.exists('../html/major_keywords'):
    os.mkdir('../html/major_keywords')

if not os.path.exists('../html/map'):
    os.mkdir('../html/map')

if not os.path.exists('../html/country'):
    os.mkdir('../html/country')

if not os.path.exists('../html/city'):
    os.mkdir('../html/city')

if not os.path.exists('../html/metrics'):
    os.mkdir('../html/metrics')

if not os.path.exists('../html/wordcloud'):
    os.mkdir('../html/wordcloud')

if not os.path.exists('../html/abstractwordcloud'):
    os.mkdir('../html/abstractwordcloud')

if not os.path.exists('../html/authornetwork'):
    os.mkdir('../html/authornetwork')

if not os.path.exists('../html/errorlog'):
    os.mkdir('../html/errorlog')


# Set up the logging. Level can be DEBUG|.....
logging.basicConfig(filename='../logs/papers.log',
                    filemode='w',
                    level=logging.DEBUG)


###########################################################
# Get the papers. This will get all the metadata and store
# it in a cache directory.
# papers will be the giant LIST that has all the papers in it, each as a dictionary
papers = []

# commenting out the get stuff as my assumption is that hughs work
# will join this up
# get.get.get(pmids, papers)

# Get list of files in merged directory
merged_files_list = listdir('../cache/processed/merged/')
merged_files_list.sort()
# merged_files_list = merged_files_list[0:10]
print str(len(merged_files_list))+' merged papers to load'

# Open each one and add to papers object
for this_merged_file in merged_files_list:
    with open('../cache/processed/merged/'+this_merged_file) as fo:
        # Will be a dictionary
        this_paper = json.load(fo)
        this_paper['filename'] = this_merged_file
        papers.append(this_paper)

# input_file = 'sample_data_object'
# input_file = 'data-alspac-all-pubmed-merged-format'
# with open('../cache/raw/'+input_file) as fo:
#     papers = json.load(fo)


print str(len(papers))+' papers to process'

# exit()

###########################################################
# Clean the data - e.g. tidy institute names
clean.clean.pre_clean(papers)
clean.clean.clean_institution(papers)
# clean.clean.do_deltas(papers)

###########################################################
# Add some extra data in - i.e. geocodes and citations
add.geocode.geocode(papers, error_log)

if update_citations:
    add.citations.citations(papers, scopus_api_key, scopus_citation_max_life)

file_name = '../data/summary_added_to'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()


bibliography.bibtex.bibtex(papers)

###########################################################
# Do some actual analysis on the data. This will result in
# some CSV type files that can be analysed.
analyse.analysis.journals(papers)

# pp.pprint(papers)

analyse.analysis.abstracts(papers)
network = analyse.analysis.authors(papers)
analyse.analysis.first_authors(papers)
analyse.analysis.inst(papers)
analyse.analysis.mesh(papers)
analyse.analysis.output_csv(papers)


###########################################################
# Make some web pages

cohort_rating = html.build_htmlv2.build_home(papers)
html.build_htmlv2.build_papers(papers)
html.build_htmlv2.build_mesh(papers)
html.build_htmlv2.build_google_map(papers)
html.build_htmlv2.build_country_map(papers)
html.build_htmlv2.build_metrics(papers, cohort_rating)
html.build_htmlv2.build_abstract_word_cloud(papers)
html.build_htmlv2.build_author_network(papers, network)
html.build_htmlv2.build_error_log(papers, error_log)

# Time Log
end_time = time.time()
elapsed_time = end_time - start_time
print "End Time: " + str(end_time)
print "Elapsed Time: " + str(elapsed_time)
