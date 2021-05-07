#!/usr/bin/env python3

# print(__file__)
# print(__name__)
# print(__module__)

import sys
import os
import re
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Bio import Entrez
from pprint import pprint
import config.config as config
from pyzotero import zotero

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

    # Connect to zotero and see what's on it
    zot = zotero.Zotero(config.zotero_id, config.zotero_type, config.zotero_api_key,)

    print("Items in zotero: " + str(zot.count_items()))

    #zotero_data = zot.items(limit=10)
    zotero_data = zot.everything(zot.items())

    # pprint(zotero_data[0])

    print("Items downloaded from zotero: " + str(len(zotero_data)))

    number_of_items = len(zotero_data)
    n = 0

    for this_item in zotero_data:

        n = n + 1
        print("\n" + str(n) + "/" + str(number_of_items))
        print(this_item['key'])
        pmid = None
        # Not interested in attachment or note itemTypes
        if this_item['data']['itemType'] in ('attachment', 'note', 'book', 'bookSection'):
            print("itemType: " + this_item['data']['itemType'] + " skipping")
            continue

        # If there is already a pmid in the extra field then skip this one
        if re.search(r'PMID:', this_item['data']['extra']):
            print("PMID already exists for " + this_item['key'])
            continue

        try:
            print("This DOI: " + this_item['data']['DOI'])

            ######################
            # search pubmed for doi and get pmid
            try:
                Entrez.email = config.pubmed_email       # Always tell NCBI who you are
                handle = Entrez.esearch(db="pubmed", term=this_item['data']['DOI']+'[Location ID]')

                pmid_data = {}
                pmid_data = Entrez.read(handle)
                handle.close()

                if pmid_data.get('Count') == "1":
                    pmid = pmid_data['IdList'][0]
                    print("PMID from PubMed: " + pmid)

            except ValueError as e:
                error = "Pubmed search for DOI: "+this_item['data']['DOI']+" error: "+str(e)
                print(error)
            except RuntimeError as e:
                error = "Pubmed search for DOI: "+this_item['data']['DOI']+" error: "+str(e)
                print(error)

            try:
                if pmid is not None:
                    
                    # Append the PMID to whatever is there
                    this_item['data']['extra'] = 'PMID:' + pmid + "\n" + this_item['data']['extra']

                    try:
                        # Note that this seems to require the whole item - just putting
                        # the extra field in fails
                        status = zot.update_item(this_item)

                        if status:
                            print("PMID added successfully (" + stc(pmid) + ")")
                    except:
                        print(status)
            except:
                print("No PMID found")
                pass
        except:
            print("No DOI for: ")
            print(this_item)




if __name__ == '__main__':
    try:
        main(sys.argv[1:])
    except KeyboardInterrupt:
        sys.exit(0)