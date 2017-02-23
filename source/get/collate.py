import papersZotero as pz
import papersCache as pc
import getDoi as pd
import getPubmed as pm
import papersMerge as pMerge
import hashlib
import re
import copy
import pprint
import json
import os
import config.config as config

def collate():
  #first, check if config.use_cached_merged_only is set to 1
  #if so, just load these and return
  if config.use_cached_merge_only == 1:
    #get list of currently merged papers
    merged_list = pc.getCacheList(filetype='/processed/merged')
    merged_papers = []
    print "use_cached_merge_only set to 1 so loading papers straight from processed/merged"
    for merged_paper in merged_list:
      print "Loading cached merged paper: "+merged_paper+"."
      merged_papers.append(pc.getCacheData(filetype='/processed/merged', filenames=[merged_paper])[merged_paper])

    return merged_papers

  #first, check for new papers from zotero repo
  zot = pz.zotPaper()

  #get list from cache
  zot_cache = pc.getCacheList(filetype='/raw/zotero')
  doi_cache = pc.getCacheList(filetype='/raw/doi')
  pm_cache = pc.getCacheList(filetype='/raw/pubmed')

  zot.collection = config.zotero_collection

  #get list of all keys in this zotero instance
  zot.getPapersKeys()
  
  new_keys = []
  new_papers = []

  #we may want to re-download the data from zotero
  #if config has the 'zotero_get_all' flag set to '1', make sure we get all papers not just new ones (i.e. load from cached file)
  if config.zotero_get_all == 1:
    new_keys = zot.papers_keys
  else:
    for num, paper_key in enumerate(zot.papers_keys):
      if paper_key not in zot_cache:
        new_keys.append(paper_key)
      else:
        #get the previously downloaded papers from the cache
        new_paper = pc.getCacheData(filetype='/raw/zotero', filenames=[paper_key])[paper_key]
        #check itemType - if it's 'note', we can ignore
        if new_paper['data']['itemType'] != 'note':
          new_papers.append(new_paper['data'])
  
  #get all new papers
  zot.getPapersList(key_list = new_keys)

  #cache the new papers
  for num, paper in enumerate(zot.papers):
    #cache zotero data
    pc.dumpJson(paper['key'], paper, 'raw/zotero')
    #add to new_papers for later doi/pubmed data retrieval
    #check itemType - if it's 'note', we can ignore
    if paper['data']['itemType'] != 'note':
      new_papers.append(paper['data'])

  #now check new_papers for doi or pubmed id and retrieve if required
  for paper in new_papers:
    print 'Getting doi/pubmed data for paper: '+paper['title']+' (zotero key: '+paper['key']+')'
    paper['doi_data'] = {}
    paper['pmid_data'] = {}
    if 'DOI' in paper and paper['DOI'] != "":
      #check if this is the full url, and if so, strip the parent
      #note that this is done again in getDoi.getDoi, but will not find the filename in the cache if we don't process this
      check_doi = re.match(r'^http://dx\.doi\.org/', paper['DOI'])
      if check_doi is not None:
        doi = re.sub(r'^http://dx\.doi\.org/', '', doi)
      #as doi's use '/' chars, we do an md5 of the doi as the filename
      doi_filename = hashlib.md5(paper['DOI']).hexdigest()

      #check if paper data in doi cache (only if config.use_pubmed_doi_cache is not 1
      if doi_filename not in doi_cache or config.use_doi_pubmed_cache != 1:
        doi_paper = pd.getDoi(paper['DOI'])
        paper['doi_data'] = doi_paper
        #data is automatically cached by getDoi
      else:
        paper['doi_data'] = pc.getCacheData(filetype='/raw/doi', filenames=[doi_filename])[doi_filename]

    #get pubmed data
    if 'extra' in paper:
      pmid_matches = re.search(r'PMID: ([0-9]{1,8})', paper['extra'])
      if pmid_matches is not None:
        paper['pmid'] = pmid_matches.group(1)
        #check if paper data in pm cache (only if config.use_doi_pubmed_cache is not 1)
        if paper['pmid'] not in pm_cache or config.use_doi_pubmed_cache != 1:
          pm_paper = pm.getPubmed(paper['pmid'])
          paper['pmid_data'] = pm_paper
          #data is automatically cached by getPubmed
        else:
          paper['pmid_data'] = pc.getCacheData(filetype='/raw/pubmed', filenames=[paper['pmid']])[paper['pmid']]
      
  #now do merge data
  merged_papers = {}
  #get list of currently merged papers
  merged_list = pc.getCacheList(filetype='/processed/merged')
  #merged_papers = []
  template_file = open(config.config_dir+'/data-doi-template', 'r')
  template = json.load(template_file)
  template_file.close()
  #template = pc.getCacheData(filenames=['data-doi-template'])['data-doi-template']
  for paper in new_papers:
    merged_paper = {}
    doi_data = paper['doi_data']
    pmid_data = paper['pmid_data']
    #delete these from paper
    del paper['doi_data']
    del paper['pmid_data']
    #create new filename
    filename = hashlib.md5(paper['title'].encode('ascii', 'ignore')).hexdigest() #md5 of title

    #if we aren't set to merge all papers, ignore existing files
    if config.merge_all != 1:
      if filename in merged_list:
        print "Merged file: "+filename+" already exists. Ignoring as merge_all not set to 1 in config.ini. File being loaded from cache."
        paper = pc.getCacheData(filetype='/processed/merged', filenames=[filename])[filename]
        continue

    print "Merging to filename: "+filename

    mgr = pMerge.Merge()
    map_file = open(config.config_dir+'/data-pubmed-doi-jsonpath', 'r')
    mgr.mapping = json.load(map_file)
    map_file.close()

    #mgr.mapping = pc.getCacheData(filenames=['data-pubmed-doi-jsonpath'])['data-pubmed-doi-jsonpath']
    mgr.src = pmid_data
    mgr.dest = copy.deepcopy(template)
    mgr.mapSrc()
    merged_paper['pmid_data'] = copy.deepcopy(mgr.dest)

    #now do zotero data
    map_file = open(config.config_dir+'/data-zotero-doi-jsonpath', 'r')
    mgr.mapping = json.load(map_file)
    map_file.close()
    #mgr.mapping = pc.getCacheData(filenames=['data-zotero-doi-jsonpath'])['data-zotero-doi-jsonpath']
    mgr.src = paper
    mgr.dest = copy.deepcopy(template)
    mgr.mapSrc()
    merged_paper['zotero_data'] = copy.deepcopy(mgr.dest)

    #now merge doi_data, merged_paper
    #starting with zotero_data, we merge using an empty mapping and setting the dest as zotero_data and src as pmid_data
    mgr.mapping = {}
    mgr.src = merged_paper['pmid_data']
    mgr.dest = merged_paper['zotero_data']
    mgr.mapSrc()
    merged_paper['pmid_zotero_data'] = copy.deepcopy(mgr.dest)

    #now set the src to be doi_data and merge to the template (otherwise the output data gets oddly formatted)
    mgr.mapping = {}
    mgr.src = doi_data
    mgr.dest = copy.deepcopy(template)
    mgr.mapSrc()
    merged_paper['doi_data'] = copy.deepcopy(mgr.dest)

    #now set the src to be doi_data and merge to the dest (pmid/zotero data)
    mgr.mapping = {"$.DOI": "$.DOI", "$.title": "$.title"}
    mgr.dest = merged_paper['doi_data']
    mgr.src = merged_paper['pmid_zotero_data']
    mgr.mapSrc()

    merged_papers[filename] = copy.deepcopy(mgr.dest)
    #merged_papers.append(copy.deepcopy(mgr.dest))
    pc.dumpJson(filename, merged_papers[filename], 'processed/merged')

  return merged_papers 

if __name__ == "__main__":
  print "Collate data"
  
  collate()

