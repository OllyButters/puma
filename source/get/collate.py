import papersZotero as pz
import papersCache as pc
import getDoi as pd
import getPubmed as pm
import papersMerge as pMerge
import hashlib
import re
import copy
import pprint

def main():
  #first, check for new papers from zotero repo
  zot = pz.zotPaper()

  #get list from cache
  zot_cache = pc.getCacheList(filetype='zotero')
  doi_cache = pc.getCacheList(filetype='doi')
  pm_cache = pc.getCacheList(filetype='pubmed')

  zot.collection = 'ALSPAC_PAPERS_ALL'

  zot.getPapersList()

  new_papers = []

  for num, paper in enumerate(zot.papers):
    if paper['key'] not in zot_cache:
      #cache zotero data
      pc.dumpJson(paper['key'], paper, 'zotero')
      #add to new_papers for later doi/pubmed data retrieval
      new_papers.append(paper['data'])
    else:
      #just for testing ONLY otherwise we'll end up getting new data all the time
      new_paper = pc.getCacheData(filetype='zotero', filenames=[paper['key']])[paper['key']]
      #check itemType - if it's 'note', we can ignore
      if new_paper['data']['itemType'] != 'note':
        new_papers.append(new_paper['data'])

  #now check new_papers for doi or pubmed id and retrieve if required
  for paper in new_papers:
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

      #check if paper data in doi cache
      if doi_filename not in doi_cache:
        doi_paper = pd.getDoi(paper['DOI'])
        paper['doi_data'] = doi_paper
        #data is automatically cached by getDoi
      else:
        paper['doi_data'] = pc.getCacheData(filetype='doi', filenames=[doi_filename])[doi_filename]

    #get pubmed data
    if 'extra' in paper:
      pmid_matches = re.search(r'PMID: ([0-9]{7})', paper['extra'])
      if pmid_matches is not None:
        paper['pmid'] = pmid_matches.group(1)
        #check if paper data in pm cache
        if paper['pmid'] not in pm_cache:
          pm_paper = pm.getPubmed(paper['pmid'])
          paper['pmid_data'] = pm_paper
          #data is automatically cached by getPubmed
        else:
          paper['pmid_data'] = pc.getCacheData(filetype='pubmed', filenames=[paper['pmid']])[paper['pmid']]
      
  #now do merge data
  merged_papers = {}
  #merged_papers = []
  template = pc.getCacheData(filenames=['data-doi-template'])['data-doi-template']
  for paper in new_papers:
    merged_paper = {}
    doi_data = copy.deepcopy(paper['doi_data'])
    pmid_data = copy.deepcopy(paper['pmid_data'])
    #delete these from paper
    del paper['doi_data']
    del paper['pmid_data']
    #create new filename
    filename = hashlib.md5(paper['title'].encode('ascii', 'ignore')).hexdigest() #md5 of title

    mgr = pMerge.Merge()
    mgr.mapping = pc.getCacheData(filenames=['data-pubmed-doi-jsonpath'])['data-pubmed-doi-jsonpath']
    mgr.src = pmid_data
    mgr.dest = copy.deepcopy(template)
    mgr.mapSrc()
    merged_paper['pmid_data'] = copy.deepcopy(mgr.dest)

    #now do zotero data
    mgr.mapping = pc.getCacheData(filenames=['data-zotero-doi-jsonpath'])['data-zotero-doi-jsonpath']
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
    mgr.src = doi_data
    mgr.dest = copy.deepcopy(template)
    mgr.mapSrc()
    merged_paper['doi_data'] = copy.deepcopy(mgr.dest)

    #now set the src to be doi_data and merge to the dest (pmid/zotero data)
    mgr.src = merged_paper['doi_data']
    mgr.dest = merged_paper['pmid_zotero_data']
    mgr.mapSrc()

    merged_papers[filename] = copy.deepcopy(mgr.dest)
    #merged_papers.append(copy.deepcopy(mgr.dest))
    pc.dumpJson(filename, merged_papers[filename], 'merged')

  return merged_papers 

