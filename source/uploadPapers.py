import get.papersCache as pc
import get.papersZotero as pz
import get.getDoi as getDoi
import sys
import getopt
import config.config as config
from Bio import Entrez
from pprint import pprint
import os
from setup.setup import build_file_tree
from shutil import copy

# upload papers to a zotero lib
# options:
# --input the input file - a json file containing a list of objects,
# each having a DOI - e.g. [{"DOI": "10.12123/hdjs.103"},{"DOI":  ... }, ...]
# --type currently only supports 'doi'
# a config.ini file in config /must/ be present with the correct Zotero api
# key and library identifier
#
# the script will setup the folder structure and then copy the input file to
# the upload folder in cache
# it will then attempt to download the doi metadata for each paper, cache
# it and then upload to Zotero
# logs are stored in [cache_dir]/processed/upload
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

    # setup dir structure if not already present
    build_file_tree()
    if not os.path.exists(os.path.join(config.cache_dir, 'upload')):
        os.mkdir(os.path.join(config.cache_dir, 'upload'))
    if not os.path.exists(os.path.join(config.cache_dir, 'processed', 'upload')):
        os.mkdir(os.path.join(config.cache_dir, 'processed', 'upload'))

    print opts

    for opt, arg in opts:
        if opt in ('i', '--input'):
            cache_file_path = arg
        elif opt in ('t', '--type'):
            src_type = arg
        #elif opt in ('c', '--collection'):
        #    collection = arg

    if cache_file_path is not None and src_type is not None:
        # copy cache file to upload dir and read in
        cache_file = os.path.split(cache_file_path)[-1]
        if os.path.isfile(cache_file_path):
            copy(cache_file_path, os.path.join(config.cache_dir, 'upload', cache_file))
        else:
            print 'Cache file does not exist'
            sys.exit(2)
        papers = pc.getCacheData(filetype='upload', filenames = [cache_file,])[cache_file]
        print papers
    else:
        sys.exit(2)

    zot = pz.zotPaper()
    #zot.collection = collection

    zot_papers = []

    get_errors = []
    #zot_upload_errors = []

    n = 1

    for paper in papers:
        print "Processing paper "+str(n)
        extra_ids = {}
        if 'Notes' in paper:
            paper['notes'] = paper['Notes']
            extra_ids = {val.split(':')[0]: val.split(':')[1] for val in paper['notes'].split('\n') if len(val.split(':')) > 1}

        #get doi data if relevant
        if 'DOI' in paper:
            print 'Found DOI: '+paper['DOI']
            print 'Getting data...'
            paper_doi = getDoi.getDoi(paper['DOI'])
            if paper_doi is not None:
                print 'Data retrieved for DOI: '+paper['DOI']
                #print paper
                paper = paper_doi
            else:
                error = "DOI: "+paper['DOI']+" data not obtained"
                get_errors.append(error)
                print error

            #search pubmed for doi and get pmid
            try:
                Entrez.email = config.pubmed_email       # Always tell NCBI who you are
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
