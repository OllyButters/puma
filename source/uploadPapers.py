import get.papersCache as pc
import get.papersZotero as pz
import get.getDoi as getDoi
import sys
import getopt
import config.config as config
from Bio import Entrez
import os

def main(argv):
  # Lets figure out some paths that everything is relative to
  # global root_dir
  path_to_papers_py = os.path.abspath(sys.argv[0])
  root_dir = os.path.dirname(os.path.dirname(path_to_papers_py))
  print 'Root directory = ' + root_dir
  config.build_config_variables(root_dir)

  try:
    opts, args = getopt.getopt(argv, "ifc", ["input=", "type=", "collection="])
  except:
    sys.exit(2)

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
  zot.collection = collection

  zot_papers = []

  n = 1

  for paper in papers:
    print "Processing paper "+str(n)
    notes = ''
    
    #get doi data if relevant
    if 'DOI' in paper:
      paper_doi = getDoi.getDoi(paper['DOI'])
      if paper_doi is not None:
        paper = paper_doi

      print paper
      #search pubmed for doi and get pmid
      try:
        Entrez.email = config.pubmed_email     # Always tell NCBI who you are
        handle = Entrez.esearch(db="pubmed", term=paper['DOI'])

        pmid_data = Entrez.read(handle)

        if pmid_data.get('Count') == "1":
          pmid = pmid_data['IdList'][0]
          notes += 'PMID: '+pmid+';'
      except ValueError, e:
        print "Pubmed search for DOI: "+paper['DOI']+" error: "+str(e)
      except RuntimeError, e:
        print "Pubmed search for DOI: "+paper['DOI']+" error: "+str(e)

    zot_papers.append(zot.mapFields(paper, src_type = src_type))
    zot_papers[-1]['notes'] = notes
    pc.dumpFile(filename=cache_file+'-'+str(n)+'-zotero-upload', data=zot_papers, filetype='processed/upload')
    print zot_papers[-1].get('title')
    n += 1

  pc.dumpFile(filename=cache_file+'-zotero-upload', data=zot_papers, filetype='processed')

  max_upload = 10

  zot_papers_chunked = [zot_papers[i:i+max_upload] for i in range(0, len(zot_papers), max_upload)]
  zot_papers_uploaded = []

  for paper_chunk in zot_papers_chunked:
    try:
      print "Uploading "+str(len(paper_chunk))+" papers to collection "+zot.collection
      created = zot.create_items(paper_chunk)
      zot_papers_uploaded += created
    except Exception, e:
      print str(e)
      sys.exit(2)
  else:
    print "Success\n"
    print str(len(zot_papers_uploaded))+" papers added"
    

if __name__ == '__main__':
  main(sys.argv[1:])
