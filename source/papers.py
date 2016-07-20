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
import html.build_htmlv2
import bibliography.bibtex

##########################################################
# Get all the paper metadata from pubmed and do stuff with it
##########################################################
# Starting to hack around with using generic template and not using pubmed

__author__ = "Olly Butters, Hugh Garner"
__date__ = 8/7/16
__version__ = '0.2.4'


# Options - these should get moved out into a config file
# Stick a flag into see if we want to update the citations
update_citations = True
scopus_api_key = '8024d746590aade6be6856a22a734783'
scopus_citation_max_life = 30  # days


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

if (os.path.exists('../cache/geodata') is False):
    os.mkdir('../cache/geodata')

if (os.path.exists('../data') is False):
    os.mkdir('../data')

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

if not os.path.exists('../html/metrics'):
    os.mkdir('../html/metrics')

if not os.path.exists('../html/wordcloud'):
    os.mkdir('../html/wordcloud')


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

# input_file = 'sample_data_object'
input_file = 'data-alspac-all-pubmed-merged-format'
with open('../cache/raw/'+input_file) as fo:
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

# analyse.analysis.abstracts(pmids, papers)
analyse.analysis.authors(papers)
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
html.build_htmlv2.build_metrics(papers, cohort_rating)

# html.build_html.build_yearly(papers)
# html.build_html.build_mesh(papers)
# html.build_html.build_summary(papers)
# html.build_html.build_google_map(papers)
html.build_html.build_google_heat_map()
