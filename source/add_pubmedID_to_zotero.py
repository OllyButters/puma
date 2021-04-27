#!/usr/bin/env python3


import get.papersZotero as pz
import sys
import config.config as config
from Bio import Entrez
from pprint import pprint
import os


# Update zotero with a PMID where one can be got from PubMed using the PubMed API.
# This connects using the details in the config.ini file.
# This requires READ/WRITE access to the zotero library.
def main(argv):
    # Lets figure out some paths that everything is relative to
    # global root_dir
    path_to_papers_py = os.path.abspath(sys.argv[0])
    root_dir = os.path.dirname(os.path.dirname(path_to_papers_py))
    print('Root directory = ' + root_dir)

    config.build_config_variables(root_dir)

    # instantiate the class
    zot = pz.zotPaper()

    zot.collection = config.zotero_collection

    # get list of ALL keys in this zotero library instance
    zot.getPapersKeys()

    new_keys = []
    new_keys = zot.papers_keys

    # get all new papers. Note this will not write any to disk until it has got all of them.
    zot.getPapersList(key_list=new_keys)

    for num, paper in enumerate(zot.papers):
        print(num)
        pprint(paper)

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
                print("PMID data: " + pmid)

        except ValueError as e:
            error = "Pubmed search for DOI: "+paper['data']['DOI']+" error: "+str(e)
            print(error)
        except RuntimeError as e:
            error = "Pubmed search for DOI: "+paper['data']['DOI']+" error: "+str(e)
            print(error)

            # Extra
            # PMID: ####

        if pmid:
            update = {}
            update['key'] = paper['key']
            update['extra'] = 'PMID:'+pmid

            zot.uploadExtra(update)


if __name__ == '__main__':
    main(sys.argv[1:])
