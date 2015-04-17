#! /usr/bin/env python

##########################################################
#Get all the paper metadata from pubmed and do stuff with it
##########################################################

#7/11/14

import csv
import json
import os.path
#import gdata.docs.service
import logging


from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are

#Look at the PubModel. See http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
override_pubmodel=False


def get(pmids, papers):
    
    ###########################################################
    #Read in the list of PMIDs from an external csv file in this directory.
    #The file must be just pmids - one per line, and no extra lines at the end.
    #pmids = []
    with open('../inputs/pmids.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            pmids.append(row[0])

        num_pmids=len(pmids)
        print str(num_pmids)+' PMIDs to process'
            
            
    ###########################################################
    #Try to build a papers dictionary with PMID as the index
        
    for this_pmid in pmids:
                
        print 'Working on '+this_pmid
        logging.info('Working on %s',this_pmid)
                
        #Build a cache of all the pmid data so we don't keep downloading it.
        #Check that cache first when looking for a PMID, if it's not there then
        #go and get the data.
        file_name='../cache/'+this_pmid 
        if not os.path.isfile(file_name):
        #This PMID data not in cache, so download it.
            logging.info('Downloading %s', this_pmid)
            handle = Entrez.efetch(db="pubmed", id=this_pmid, retmode="xml")
            #record = Entrez.read(handle)
           
            #Should check to see if anything sensible was returned - eg check size

            file_name='../cache/'+this_pmid
            fo = open(file_name, 'w')
            fo.write(handle.read())
            #fo.write(record[0])
            fo.close()
            
        #Open and parse cached file
        fo = open(file_name, 'r')
        #record = json.loads(fo.read())
        record_a = Entrez.read(fo)
        fo.close()
        
        #Should only be one record in each, so just grab that.
        record = record_a[0]
        #print record

                    
        #Define the info we want for this paper
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
            this_paper['MeshHeadingList'] = []
            for this_mesh_heading in record['MedlineCitation']['MeshHeadingList']:
                temp = {'DescriptorName':this_mesh_heading['DescriptorName'], 'MajorTopicYN':this_mesh_heading['DescriptorName'].attributes['MajorTopicYN']}
                this_paper['MeshHeadingList'].append(temp)
        except:
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

        #Look at the PubModel. See http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
        try:
            this_paper['PubModel'] = record['MedlineCitation']['Article'].attributes['PubModel']
        except:
            print 'No PubModel listed!'
            logging.info('No Pubmodel listed')
            logging.warn('No Pubmodel text')

        #Get the year based on the pubmodel
        #Look at the PubModel. See http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
        #Can decide if this is respected with override_pubmodel
        if this_paper['PubModel'] == 'Print-Electronic' and override_pubmodel:
            try:
                this_paper['Year'] = record['MedlineCitation']['Article']['ArticleDate'][0]['Year']
            except:
                print 'No ArticleDate Year listed!'
                logging.info('No ArticleDate Year listed')
                logging.warn('No ArticleDate Year text')
        else:
            #Publication date Year
            try:
                this_paper['Year'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
            except:
                try:
                    #Could be that it is a date range see http://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#articledate
                    temp = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate']
                    #lets assume that the first 4 characters are the year....
                    this_paper['Year']=temp[0:4]
                except:
                    print 'No PubDate Year listed!'
                    logging.info('No PubDate Year listed')
                    logging.warn('No PubDate Year text')

        #Try the same as above for the month
        if this_paper['PubModel'] == 'Print-Electronic' and override_pubmodel:
            try:
                this_paper['Month'] = record['MedlineCitation']['Article']['ArticleDate'][0]['Month']
            except:
                print 'No ArticleDate Month listed!'
                logging.info('No ArticleDate Month listed')
                logging.warn('No ArticleDate Month text')
        else:
            #Publication date Year
            try:
                this_paper['Month'] = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
            except:
                print 'No PubDate Month listed!'
                logging.info('No PubDate Month listed')
                logging.warn('No PubDate Month text')

    #Publication date Month
#    try:
#        this_paper['ArticleDateMonth'] = record['MedlineCitation']['Article']['ArticleDate'][0]['Month']
#    except:
#        print 'No ArticleDateMonth listed!'
#        logging.info('No ArticleDateMonth listed')
#        logging.warn('No ArticleDateMonth text')

        #Add an extra field that we will fill up with other data
        this_paper['Extras'] = {}


        #Add this_paper info to the main dict
        papers[this_pmid]=this_paper


    #Save it for later
    file_name='../data/summary'
    fo = open(file_name, 'wb')
    fo.write(json.dumps(papers, indent=4))
    fo.close()

#print '####All papers'
#print papers
