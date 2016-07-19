import csv
import json
import os.path
import logging
import papersCache as pc


from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are

#Look at the PubModel. See http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
override_pubmodel=False

def getPubmed(this_pmid):
  print 'Working on '+this_pmid
  logging.info('Working on %s',this_pmid)

  logging.info('Downloading %s', this_pmid)
  handle = Entrez.efetch(db="pubmed", id=this_pmid, retmode="xml")

  #pmid_xml_data = handle.read()

  #pc.dumpFile(this_pmid+'.xml', pmid_xml_data, 'pubmed/raw')

  #pmid_data = Entrez.read(pmid_xml_data)[0]
  try:
    pmid_data = Entrez.read(handle)[0]
    pc.dumpJson(this_pmid, pmid_data, 'pubmed')
    return pmid_data
  except:
    logging.warn('Unable to read pmid %s', this_pmid)
    return None
  

