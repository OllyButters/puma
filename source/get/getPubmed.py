import csv
import json
import os.path
import logging
import papersCache as pc
from pprint import pprint

from Bio import Entrez
Entrez.email = "olly.butters@bristol.ac.uk"     # Always tell NCBI who you are

#Look at the PubModel. See http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
override_pubmodel=False

def getPubmed(this_pmid):
  print 'Working on '+this_pmid
  logging.info('Working on %s',this_pmid)

  logging.info('Downloading %s', this_pmid)
  handle = Entrez.efetch(db="pubmed", id=this_pmid, retmode="xml")

  pmid_xml_data = handle.read()

  xml_file_loc = pc.dumpFile(this_pmid+'.xml', pmid_xml_data, 'raw/pubmed/xml')


  try:
    xml_file = open(xml_file_loc, 'r')
    pmid_data = Entrez.read(xml_file)[0]
    xml_file.close()

    ###
    # data processing
    # some dp is required as Entrez.read returns a subclassed Dict type with additonal xml data as attributes. These are not serialised by the json dump so we need to include them in another way.
    ###

    #add asterisk to major mesh headings
    for this_mesh_heading in pmid_data['MedlineCitation']['MeshHeadingList']:
      this_mesh_heading['MajorTopicYN'] = this_mesh_heading['DescriptorName'].attributes['MajorTopicYN']

    ###
    #
    # end data processing
    #
    ###

    pc.dumpJson(this_pmid, pmid_data, 'raw/pubmed')
    return pmid_data
  except:
    logging.warn('Unable to read pmid %s', this_pmid)
    return None
  

