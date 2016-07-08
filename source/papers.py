#! /usr/bin/env python

import json
import os.path
import logging
# import pprint

# import get.get
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

__author__ = "Olly Butters, Hugh Garner"
__date__ = 8/7/16
__version__ = '0.2.4'


# Options
# Stick a flag into see if we want to update the citations
update_citations = True

# pp = pprint.PrettyPrinter(indent=4)

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

with open('../cache/raw/sample_data_object') as fo:
    papers = json.load(fo)

print str(len(papers))+' papers to process'


###########################################################
# Clean the data - e.g. tidy institute names
clean.clean.pre_clean(papers)
clean.clean.clean_institution(papers)
# clean.clean.do_deltas(papers)

###########################################################
# Add some extra data in - i.e. geocodes and citations
add.geocode.geocode(papers)

if update_citations:
    add.citations.citations(papers)

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

# analyse.analysis.abstracts(pmids, papers)
analyse.analysis.authors(papers)
analyse.analysis.first_authors(papers)
analyse.analysis.inst(papers)
analyse.analysis.mesh(papers)
analyse.analysis.output_csv(papers)


###########################################################
# Make some web pages
html.build_html.build_yearly(papers)
html.build_html.build_mesh(papers)
html.build_html.build_summary(papers)
html.build_html.build_google_map(papers)
html.build_html.build_google_heat_map()
