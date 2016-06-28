#! /usr/bin/env python

##########################################################
#Get all the paper metadata from pubmed and do stuff with it
##########################################################
#Starting to hack around with using generic template and not using pubmed

#need to have a well defined template of the data model to refer to.

__author__ = "Olly Butters, Hugh Garner"
__date__ = 28/6/16
__version__ = '0.2.1'

import csv
import json
import os.path
#import gdata.docs.service
import logging

import get.get
import clean.clean
import add.geocode
import add.citations
import analyse.analysis
import html.build_html
import bibliography.bibtex

#Stick a flag into see if we want to update the citations
update_citations = True


###########################################################
#Make sure the directory structure is set up first.
#Everything in the cache is grabbed from elsewhere, or built on the fly,
#so it should all be considerd ready to be deleted at any point!
if (os.path.exists('../cache') == False):
    os.mkdir('../cache')

#The raw, unprocessed data.
if (os.path.exists('../cache/raw') == False):
    os.mkdir('../cache/raw')

if (os.path.exists('../cache/raw/pubmed') == False):
    os.mkdir('../cache/raw/pubmed')

if (os.path.exists('../cache/processed') == False):
    os.mkdir('../cache/processed')

if (os.path.exists('../data') == False):
    os.mkdir('../data')

if (os.path.exists('../html') == False):
    os.mkdir('../html')

if (os.path.exists('../html/mesh') == False):
    os.mkdir('../html/mesh')


#Set up the logging
logging.basicConfig(filename='../data/papers.log',filemode='w',level=logging.DEBUG)


###########################################################
#Get the papers. This will get all the metadata and store
#it in a cache directory.
#papers will be the giant object that has all the papers in it
#pmids is the list of papers
papers = {}
pmids = []

#commenting out the get stuff as my assumption is that hughs work will join this up
#get.get.get(pmids, papers)

#use a md5 hash of article title for the id. here are two taken from hughs eg
#might want to think about how we generate these hashes - should we process the
#titles a bit, e.g. get rid of punctuation that might make them a little
#ambiguous?
paper_list = ['e2cdfaede7d8e4207820a6ea36e6e01b', '47d268cdce86aa9248ea534ea6b5b5eb']
print paper_list
#exit()

#sort the list - have a better idea of the order things will run
#paper_lis = paper_list.sort()
print paper_list

#file_name='../cache/raw/'+paper_list[1]
#print file_name

#with open(file_name) as fo:
#    papers=json.load(fo)

#print papers
#print '============='
#print papers[0]['author'][0]['affiliation'][0]['name']
#exit()

#file_name='../data/summary'
#with open(file_name) as fo:
#    papers=json.load(fo)

print str(len(paper_list))+' papers to process'
#print str(len(papers))+' papers processed'

###########################################################
#Clean the data - e.g. tidy institute names
clean.clean.pre_clean(paper_list)
clean.clean.clean_institution(paper_list)
#clean.clean.do_deltas(papers)

#exit()

#Save it for later
#file_name='../data/summary_cleaned'
#fo = open(file_name, 'wb')
#fo.write(json.dumps(papers, indent=4))
#fo.close()

###########################################################
#Add some extra data in - i.e. geocodes and citations
add.geocode.geocode(paper_list)

exit()

if update_citations:
    add.citations.citations(pmids,papers)


file_name='../data/summary_added_to'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()


bibliography.bibtex.bibtex(pmids,papers)

###########################################################
#Do some actual analysis on the data. This will result in
#some CSV type files that can be analysed.
analyse.analysis.journals(pmids, papers)
#analyse.analysis.abstracts(pmids, papers)
analyse.analysis.authors(pmids, papers)
analyse.analysis.first_authors(pmids, papers)
analyse.analysis.inst(pmids, papers)
analyse.analysis.mesh(pmids, papers)

analyse.analysis.output_csv(pmids, papers)

###########################################################
#Make some web pages
html.build_html.build_yearly(pmids, papers)
html.build_html.build_mesh(pmids, papers)
html.build_html.build_summary(pmids, papers)
html.build_html.build_google_map(pmids, papers)
html.build_html.build_google_heat_map()
