#! /usr/bin/env python

##########################################################
#Get all the paper metadata from pubmed and do stuff with it
##########################################################

#21/9/14

import csv
import json
import os.path
#import gdata.docs.service

import tools
import analysis

from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are


###########################################################
#Make sure the directory structure is set up first
if (os.path.exists('../cache') == False):
    os.mkdir('../cache')

if (os.path.exists('../data') == False):
    os.mkdir('../data')




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


    #print all
    #print '#####ALL'
    #print record


    #print '#####Authors'
    #print all authors
    #print record[0]['MedlineCitation']['Article']['AuthorList']
    
    
    #print out lastname and affiliation for each author
    #for this_one in record[0]['MedlineCitation']['Article']['AuthorList']:
     #   print this_one['LastName']
      #  print this_one['Affiliation']
               
    #all the info for this paper
    this_paper = {}
    this_paper['pmid']=this_pmid
    this_paper['ArticleTitle'] = record['MedlineCitation']['Article']['ArticleTitle']
    this_paper['Journal'] = record['MedlineCitation']['Article']['Journal']['ISOAbbreviation']

    try:
        this_paper['AbstractText'] = record['MedlineCitation']['Article']['Abstract']['AbstractText']
    except:
        print 'No Abstract text'
    
    #Mesh keywords
    try:
        this_paper['MeshHeadingList'] = record['MedlineCitation']['MeshHeadingList']
    except:
        #Should log this
        print 'No mesh headings'

    #Author list
    try:
        this_paper['AuthorList'] = record['MedlineCitation']['Article']['AuthorList']
    except:
        print 'No Authors listed!'

    
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
