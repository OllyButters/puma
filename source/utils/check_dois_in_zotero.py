#!/usr/bin/env python3


# import get.papersZotero as pz
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config.config as config
# from Bio import Entrez
# from pprint import pprint

from pyzotero import zotero
import csv


# run like:
# ./check_dois_in_zotero.py --config config_whatever.ini
# This connects using the details in the config.ini file.
# All DOIs get cast to lowercase to enable matching.
def main(argv):
    # Lets figure out some paths that everything is relative to
    # global root_dir
    path_to_this_file = os.path.abspath(sys.argv[0])
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(path_to_this_file)))
    print('Root directory = ' + root_dir)

    config.build_config_variables(root_dir)

    doi_list_path = "/home/olly/phps_dois.csv"

    # Get the list of DOIs in the input file
    input_dois = []
    with open(doi_list_path) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            input_dois.append(row[0].strip().lower())
            # print(row)

    print(str(len(input_dois)) + " DOIs in input file")
#    doi_list_duplicates = input_dois - list(set(input_dois))
    #doi_list_duplicates = list(input_dois.difference(set(input_dois)))
    # doi_list_duplicates = list(set([x for x in input_dois if input_dois.count(x) > 1]))
    doi_list_duplicates = [item for item in set(input_dois) if input_dois.count(item) > 1]
    print(str(len(doi_list_duplicates)) + " DUPLICATE DOIs in input file")
    print(str(len(set(input_dois))) + " UNIQUE DOIs in input file. \n Duplicates:")
    print(doi_list_duplicates)


    # Connect to zotero and see what's on it
    zot = zotero.Zotero(config.zotero_id, config.zotero_type, config.zotero_api_key,)

    print("Items in zotero: " + str(zot.count_items()))

    # zotero_data = zot.items()
    zotero_data = zot.everything(zot.items())

    print("Items downloaded from zotero: " + str(len(zotero_data)))

    existing_dois = []
    ignored_itemtypes = 0
    for this_item in zotero_data:
        try:
            # print(this_item['data']['DOI'])
            existing_dois.append(this_item['data']['DOI'].lower())
        except:
            # If there are itemType = book, bookSection then the DOI
            # gets put in the Extra field. PUMA ignores these anyway, so 
            # don't worry too much about them here either
            if this_item['data']['itemType'] in ('book', 'bookSection'):
                ignored_itemtypes = ignored_itemtypes + 1

                # Pretty poor way of adding book DOIs to list, assumes
                # DOI is the only item in the extra field. 
                this_extra = this_item['data']['extra']
                if "DOI:" in this_extra:
                    this_extra = this_extra.replace("DOI:", "").strip()
                    existing_dois.append(this_extra.lower())


                continue
            print(this_item)

    print("Items downloaded from zotero with DOIs (including books): " + str(len(existing_dois)))

    doi_diff = list(set(input_dois) - set(existing_dois))
    print(str(len(doi_diff)) + " missing DOIs")
    # print(str(ignored_itemtypes) + " books/bookSection in missing DOIs")

    print("Missing DOIs:")
    for this_item in sorted(doi_diff):
        print(this_item)


    doi_diff = list(set(existing_dois) - set(input_dois))
    print(str(len(doi_diff)) + " missing DOIs")
    # print(str(ignored_itemtypes) + " books/bookSection in missing DOIs")

    print("Missing DOIs:")
    for this_item in sorted(doi_diff):
        print(this_item)


if __name__ == '__main__':
    main(sys.argv[1:])
