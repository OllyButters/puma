import csv
import json
import os.path
import logging
import papersCache as pc


from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are

#Look at the PubModel. See http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
override_pubmodel=False

def getPubmed(pmid):
  print 'Working on '+this_pmid
  logging.info('Working on %s',this_pmid)

  logging.info('Downloading %s', this_pmid)
  handle = Entrez.efetch(db="pubmed", id=this_pmid, retmode="xml")

  pmid_xml_data = handle.read()

  pc.dumpFile(pmid+'.xml', pmid_xml_data, 'pubmed/raw')

  pmid_data = Entrez.read(pmid_xml_data)[0]
  
  pc.dumpJson(pmid, pmid_data, 'pubmed')

  return pmid_data
