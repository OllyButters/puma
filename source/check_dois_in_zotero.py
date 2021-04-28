#!/usr/bin/env python3


import get.papersZotero as pz
import sys
import config.config as config
from Bio import Entrez
from pprint import pprint
import os
from pyzotero import zotero
import csv

# This connects using the details in the config.ini file.
def main(argv):
    # Lets figure out some paths that everything is relative to
    # global root_dir
    path_to_papers_py = os.path.abspath(sys.argv[0])
    root_dir = os.path.dirname(os.path.dirname(path_to_papers_py))
    print('Root directory = ' + root_dir)

    config.build_config_variables(root_dir)

    doi_list_path = "/home/olly/temp.csv"


    # Get the list of DOIs in the input file
    input_dois = []
    with open(doi_list_path) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            input_dois.append(row[0])
            # print(row)

    print(str(len(input_dois)) + " DOIs in input file")

    # Connect to zotero and see what's on it
    zot = zotero.Zotero(config.zotero_id, config.zotero_type, config.zotero_api_key,)

    print("Items in zotero: " + str(zot.count_items()))

    zotero_data = zot.items()

    existing_dois = []
    for this_item in zotero_data:
        existing_dois.append(this_item['data']['DOI'])

    print("Items downloaded from zotero: " + str(len(existing_dois)))
 
    doi_diff = list(set(input_dois) - set(existing_dois))
    print(str(len(doi_diff)) + " missing DOIs")
 
    print("Missing DOIs:")
    for this_item in doi_diff:
        print(this_item)



if __name__ == '__main__':
    main(sys.argv[1:])
