#! /usr/bin/env python

##########################################################
#Get all the paper metadata from pubmed and do stuff with it
##########################################################

#17/10/14

import csv
import json
import os.path
#import gdata.docs.service
import logging


import tools
import analysis
import build_html

from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are


###########################################################
#Make sure the directory structure is set up first
if (os.path.exists('../cache') == False):
    os.mkdir('../cache')

if (os.path.exists('../data') == False):
    os.mkdir('../data')

if (os.path.exists('../html') == False):
    os.mkdir('../html')


#Set up the logging
logging.basicConfig(filename='../data/papers.log',filemode='w',level=logging.INFO)


###########################################################
#Read in the list of PMIDs from an external csv file in this directory.
#The file must be just pmids - one per line, and no extra lines at the end.
pmids = []
with open('../inputs/pmids.csv', 'rb') as csvfile:
    f = csv.reader(csvfile)
    for row in f:
        #print row[0]
        pmids.append(row[0])

#print pmids

num_pmids=len(pmids)
print str(num_pmids)+' PMIDs to process'


###########################################################
#Try to build a papers dictionary with PMID as the index
papers = {}

for this_pmid in pmids:

    print 'Working on '+this_pmid
    logging.info('Working on %s',this_pmid)

    #Build a cache of all the pmid data so we don't keep downloading it.
    #Check that cache first when looking for a PMID, if it's not there then
    #go and get the data.
    file_name='../cache/'+this_pmid 
    if not os.path.isfile(file_name):
        #This PMID data not in cache, so download it.
        handle = Entrez.efetch(db="pubmed", id=this_pmid, retmode="xml")
        record = Entrez.read(handle)
        
        #Should check to see if anything sensible was returned - eg check size

        #Save it for later
        file_name='../cache/'+this_pmid
        fo = open(file_name, 'wb')
        fo.write(json.dumps(record[0], indent=4))
        fo.close()

    #Open and parse cached file
    fo = open(file_name, 'r')
    record = json.loads(fo.read())
    fo.close()

               
    #all the info for this paper
    this_paper = {}
    this_paper['pmid']=this_pmid
    this_paper['ArticleTitle'] = record['MedlineCitation']['Article']['ArticleTitle']
    this_paper['Journal'] = record['MedlineCitation']['Article']['Journal']['ISOAbbreviation']

    #Journal volume
    try:
        this_paper['JournalVolume'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
    except:
        print 'No JournalVolume'
        logging.info('No Journal volume')
        logging.warn('No Journal volume')

    #Abstract
    try:
        this_paper['AbstractText'] = record['MedlineCitation']['Article']['Abstract']['AbstractText']
    except:
        print 'No Abstract text'
        logging.info('No abstract text')
        logging.warn('No abstract text')

    #DOI
    try:
        this_paper['doi'] = record['MedlineCitation']['Article']['ELocationID']
    except:
        print 'No doi'
        logging.info('No DOI')
        logging.warn('No DOI')
    
    #Mesh keywords
    try:
        this_paper['MeshHeadingList'] = record['MedlineCitation']['MeshHeadingList']
    except:
        #Should log this
        print 'No mesh headings'
        logging.info('No mesh headings')
        logging.warn('No abstract text')

    #Author list
    try:
        this_paper['AuthorList'] = record['MedlineCitation']['Article']['AuthorList']
    except:
        print 'No Authors listed!'
        logging.info('No Authors listed')
        logging.warn('No abstract text')

    #Publication date Year
    try:
        #this_paper['ArticleDateYear'] = record['MedlineCitation']['Article']['ArticleDate'][0]['Year']
        this_paper['Year'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
    except:
        print 'No PubDate Year listed!'
        logging.info('No PubDate Year listed')
        logging.warn('No PubDate Year text')

    #Publication date Month
#    try:
#        this_paper['ArticleDateMonth'] = record['MedlineCitation']['Article']['ArticleDate'][0]['Month']
#    except:
#        print 'No ArticleDateMonth listed!'
#        logging.info('No ArticleDateMonth listed')
#        logging.warn('No ArticleDateMonth text')


    #Add this_paper info to the main dict
    papers[this_pmid]=this_paper


#Save it for later
file_name='../data/summary'
fo = open(file_name, 'wb')
fo.write(json.dumps(papers, indent=4))
fo.close()

#print '####All papers'
#print papers

print str(num_pmids)+' PMIDs to process'
print str(len(papers))+' papers processed'


#Do some actual analysis on the data
analysis.journals(pmids, papers)
#analysis.abstracts(pmids, papers)
analysis.authors(pmids, papers)
analysis.first_authors(pmids, papers)
analysis.inst(pmids, papers)
analysis.mesh(pmids, papers)

build_html.build_yearly(pmids, papers)
build_html.build_mesh(pmids, papers)
build_html.build_summary(pmids, papers)
