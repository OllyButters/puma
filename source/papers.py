#! /usr/bin/env python

##########################################################
#Get all the paper metadata from pubmed and do stuff with it
##########################################################

#24/10/14

import csv
import json
import os.path
#import gdata.docs.service
import logging

import get.get
import clean.clean
import add.geocode
import analyse.analysis 
import html.build_html



###########################################################
#Make sure the directory structure is set up first
if (os.path.exists('../cache') == False):
    os.mkdir('../cache')

if (os.path.exists('../data') == False):
    os.mkdir('../data')

if (os.path.exists('../html') == False):
    os.mkdir('../html')

if (os.path.exists('../html/mesh') == False):
    os.mkdir('../html/mesh')


#Set up the logging
logging.basicConfig(filename='../data/papers.log',filemode='w',level=logging.WARN)


###########################################################
#Get the papers
papers = {}
pmids = []
get.get.get(pmids, papers)

print str(len(pmids))+' PMIDs to process'
print str(len(papers))+' papers processed'

###########################################################
#Clean the data - e.g. tidy institute names
clean.clean.clean_institution(pmids,papers)
clean.clean.do_deltas(papers)

#Save it for later
file_name='../data/summary_cleaned'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()



###########################################################
#Add some extra data in - e.g. geocodes 
add.geocode.geocode(pmids,papers)
file_name='../data/summary_added_to'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()


###########################################################
#Do some actual analysis on the data
analyse.analysis.journals(pmids, papers)
#analyse.analysis.abstracts(pmids, papers)
analyse.analysis.authors(pmids, papers)
analyse.analysis.first_authors(pmids, papers)
analyse.analysis.inst(pmids, papers)
analyse.analysis.mesh(pmids, papers)

###########################################################
#Make some web pages
html.build_html.build_yearly(pmids, papers)
html.build_html.build_mesh(pmids, papers)
html.build_html.build_summary(pmids, papers)
html.build_html.build_google_map(pmids, papers)
html.build_html.build_google_heat_map()
