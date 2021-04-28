#!/usr/bin/env python3

# print(__file__)
# print(__name__)
# print(__module__)

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Bio import Entrez
# from pprint import pprint
import get.papersZotero as pz
import get.papersZotero as pz
import config.config as config

# Note the sys.path stuff here is so it picks up the files in parent directories
# Need to figure out how to do this properly


# Update zotero by adding a PMID to the extra field like "pmid: 12345"
# where one can be got from PubMed using the PubMed API.
# This connects using the details in the config.ini file.
# run like:
# ./add_pubmedID_to_zotero.py --config config_whatever.ini
# This preserves what is in the extra field, only appending the pmid if
# it is missing. If it already exists then it is NOT updated.
# This requires READ/WRITE access to the zotero library.
def main(argv):
    # Lets figure out some paths that everything is relative to
    # global root_dir
    path_to_this_file = os.path.abspath(sys.argv[0])
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(path_to_this_file)))
    print('Root directory = ' + root_dir)

    config.build_config_variables(root_dir)

    # instantiate the class
    zot = pz.zotPaper()

    zot.collection = config.zotero_collection

    # get list of ALL keys in this zotero library instance
    zot.getPapersKeys()

    new_keys = []
    new_keys = zot.papers_keys

    # get all new papers. Note this will not write anything to disk
    zot.getPapersList(key_list=new_keys)

    num_papers = len(zot)
    print(str(num_papers) + " to process.")

    for num, paper in enumerate(zot.papers):
        print("\nWorking on " + str(num) + "/" + str(num_papers))
        # pprint(paper)

        # Don't care about attachment or note items
        if paper['data']['itemType'] in ('attachment', 'note'):
            print("itemType: " + paper['data']['itemType'] + " skipping")
            continue

        print(paper['key'])
        print(paper['data']['DOI'])
        print(paper['data']['extra'])

        ######################
        # search pubmed for doi and get pmid
        try:
            Entrez.email = config.pubmed_email       # Always tell NCBI who you are
            handle = Entrez.esearch(db="pubmed", term=paper['data']['DOI']+'[Location ID]')

            pmid_data = {}
            pmid_data = Entrez.read(handle)
            handle.close()

            if pmid_data.get('Count') == "1":
                pmid = pmid_data['IdList'][0]
                print("PMID from PubMed: " + pmid)

        except ValueError as e:
            error = "Pubmed search for DOI: "+paper['data']['DOI']+" error: "+str(e)
            print(error)
        except RuntimeError as e:
            error = "Pubmed search for DOI: "+paper['data']['DOI']+" error: "+str(e)
            print(error)

        try:
            pmid
            update = {}
            update['key'] = paper['key']
            update['extra'] = 'PMID:'+pmid

            zot.uploadExtra(update)
        except:
            pass


if __name__ == '__main__':
    main(sys.argv[1:])
