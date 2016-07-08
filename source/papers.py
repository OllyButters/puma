#! /usr/bin/env python

# import csv
import json
import os.path
# import gdata.docs.service
import logging

import get.get
import clean.clean
import add.geocode
import add.citations
import analyse.analysis
import html.build_html
import bibliography.bibtex

##########################################################
# Get all the paper metadata from pubmed and do stuff with it
##########################################################
# Starting to hack around with using generic template and not using pubmed

# - need to have a well defined template of the data model to refer to.
# - opening the json files in each subroutine, not that efficient.

__author__ = "Olly Butters, Hugh Garner"
__date__ = 7/7/16
__version__ = '0.2.3'


# Options
# Stick a flag into see if we want to update the citations
update_citations = True


###########################################################
# Make sure the directory structure is set up first.
# Everything in the cache is grabbed from elsewhere, or built on the fly,
# so it should all be considerd ready to be deleted at any point!
if (os.path.exists('../cache') is False):
    os.mkdir('../cache')

# The raw, unprocessed data.
if (os.path.exists('../cache/raw') is False):
    os.mkdir('../cache/raw')

if (os.path.exists('../cache/raw/pubmed') is False):
    os.mkdir('../cache/raw/pubmed')

if (os.path.exists('../cache/processed') is False):
    os.mkdir('../cache/processed')

if (os.path.exists('../data') is False):
    os.mkdir('../data')

if (os.path.exists('../html') is False):
    os.mkdir('../html')

if (os.path.exists('../html/mesh') is False):
    os.mkdir('../html/mesh')


# Set up the logging
logging.basicConfig(filename='../data/papers.log',
                    filemode='w',
                    level=logging.DEBUG)


###########################################################
# Get the papers. This will get all the metadata and store
# it in a cache directory.
# papers will be the giant object that has all the papers in it
papers = {}

# commenting out the get stuff as my assumption is that hughs work
# will join this up
# get.get.get(pmids, papers)

# use a md5 hash of article title for the id. here are two taken from hughs eg
# might want to think about how we generate these hashes - should we process
# the titles a bit, e.g. get rid of punctuation that might make them a little
# ambiguous?
paper_list = [
    'e2cdfaede7d8e4207820a6ea36e6e01b',
    '47d268cdce86aa9248ea534ea6b5b5eb']

# sort the list - have a better idea of the order things will run
paper_list.sort()
print paper_list

print str(len(paper_list))+' papers to process'

###########################################################
# Clean the data - e.g. tidy institute names
clean.clean.pre_clean(paper_list)
clean.clean.clean_institution(paper_list)
# clean.clean.do_deltas(papers)

###########################################################
# Add some extra data in - i.e. geocodes and citations
add.geocode.geocode(paper_list)

if update_citations:
    add.citations.citations(paper_list)

file_name = '../data/summary_added_to'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()


# bibliography.bibtex.bibtex(pmids,papers)

###########################################################
# Do some actual analysis on the data. This will result in
# some CSV type files that can be analysed.
analyse.analysis.journals(paper_list)

# analyse.analysis.abstracts(pmids, papers)
analyse.analysis.authors(paper_list)
analyse.analysis.first_authors(paper_list)
analyse.analysis.inst(paper_list)
analyse.analysis.mesh(paper_list)
analyse.analysis.output_csv(paper_list)

###########################################################
# Make some web pages
html.build_html.build_yearly(paper_list)
exit()
html.build_html.build_mesh(paper_list)
html.build_html.build_summary(paper_list)
html.build_html.build_google_map(paper_list)
html.build_html.build_google_heat_map()
