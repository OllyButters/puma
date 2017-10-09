import get.papersCache as pc
import get.papersZotero as pz
import get.getDoi as getDoi
import get.getPubmed as getPubmed
import copy
import json
import get.papersMerge as pMerge
import sys
import getopt
import config.config as config
from Bio import Entrez
from pprint import pprint
import os

def main(argv):
  # Lets figure out some paths that everything is relative to
  # global root_dir
  path_to_papers_py = os.path.abspath(sys.argv[0])
  root_dir = os.path.dirname(os.path.dirname(path_to_papers_py))
  print 'Root directory = ' + root_dir
  try:
    opts, args = getopt.getopt(argv, "i:t:g:", ["input=", "type=", "collection="])
  except Exception as e:
    pprint(str(e))
    sys.exit(2)
  # only pass config arg to config.ini
  sys.argv = [sys.argv[0]] + [[v, sys.argv[i+1]] for i,v in enumerate(sys.argv) if v == 'config' or v == 'c']

  config.build_config_variables(root_dir)

  pprint('Cache dir set to: '+config.cache_dir)
  cache_file = None
  src_type = None
  collection = None
  
  print opts

  for opt, arg in opts:
    if opt in ('i', '--input'):
      cache_file = arg
    elif opt in ('t', '--type'):
      src_type = arg
    elif opt in ('c', '--collection'):
      collection = arg

  if cache_file is not None and src_type is not None:
    papers = pc.getCacheData(filetype='upload', filenames = [cache_file,])[cache_file]
  else:
    sys.exit(2)

  zot = pz.zotPaper()
  #zot.collection = collection

  zot_papers = []

  get_errors = []
  zot_upload_errors = []

  n = 1

  for paper in papers:
    print "Processing paper "+str(n)
    notes = ''
    extra_ids = {}
    if 'Notes' in paper:
      paper['notes'] = paper['Notes']
      extra_ids = {val.split(':')[0]: val.split(':')[1] for val in paper['notes'].split('\n') if len(val.split(':')) > 1}
    
    #get doi data if relevant
    if 'DOI' in paper:
      paper_doi = getDoi.getDoi(paper['DOI'])
      if paper_doi is not None:
        print paper
        paper = paper_doi
      else:
        error = "DOI: "+paper['DOI']+" data not obtained"
        get_errors.append(error)
        print error 

      #search pubmed for doi and get pmid
      try:
        Entrez.email = config.pubmed_email     # Always tell NCBI who you are
        handle = Entrez.esearch(db="pubmed", term=paper['DOI']+'[Location ID]')

        pmid_data = {}

        pmid_data = Entrez.read(handle)

        handle.close()

        if pmid_data.get('Count') == "1":
          pmid = pmid_data['IdList'][0]
          try:
            paper['notes'] += 'PMID: '+pmid+';'
          except KeyError:
            paper['notes'] = 'PMID: '+pmid+';'
          extra_ids['PMID'] = pmid
      except ValueError, e:
        error = "Pubmed search for DOI: "+paper['DOI']+" error: "+str(e)
        get_errors.append(error)
        print error 
      except RuntimeError, e:
        error = "Pubmed search for DOI: "+paper['DOI']+" error: "+str(e)
        get_errors.append(error)
        print error 

    if 'PMID' in extra_ids:
      #if PMID is in extra_ids (i.e. extracted from 'Notes') get the data from pubmed. Note that this overrides the above doi data, so if both doi and pubmed present in original source data the pubmed is used; this isn't currently the case but may be in future so the pubmed data may need to be merged with the doi as per get.collate
      try:
        pubmed_data = getPubmed.getPubmed(extra_ids['PMID'])

        template_file = open(config.config_dir+'/data-doi-template', 'r')
        template = json.load(template_file)
        template_file.close()

        mgr = pMerge.Merge()
        map_file = open(config.config_dir+'/data-pubmed-doi-jsonpath', 'r')
        mgr.mapping = json.load(map_file)
        map_file.close()

        mgr.src = pubmed_data
        mgr.dest = copy.deepcopy(template)
        mgr.mapSrc()
        pubmed_paper = copy.deepcopy(mgr.dest)

        mgr.src = paper
        mgr.dest = copy.deepcopy(pubmed_paper)
        mgr.mapSrc()
        paper = copy.deepcopy(mgr.dest)

      except ValueError, e:
        error = "Pubmed search for PMID: "+extra_ids['PMID']+" error: "+str(e)
        get_errors.append(error)
        print error 
      except RuntimeError, e:
        error = "Pubmed search for PMID: "+extra_ids['PMID']+" error: "+str(e)
        get_errors.append(error)
        print error 

    zot_paper = zot.mapFields(paper, src_type = src_type)
    zot_papers.append(zot_paper)
    pc.dumpJson(filename=cache_file+'-'+str(n)+'-zotero-upload', data=zot_papers[-1], filetype='processed/upload')
    print zot_papers[-1].get('title')
    n += 1

    #dump errors to file
    pc.dumpFile(filename=cache_file+'upload-errors', filetype='processed/upload', data='\n'.join(get_errors))

    #now try uploading and dump zotero output
    try:
      print "Uploading paper "+str(n)+" to collection "+str(zot.collection)
      created = zot.create_items([zot_papers[-1]])
      pprint(created)
      pc.dumpJson(filename=cache_file+'-'+str(n)+'-zotero-upload-response', data=created, filetype='processed/upload', process=False)
    except Exception, e:
      print str(e)
      sys.exit(2)

if __name__ == '__main__':
  main(sys.argv[1:])
