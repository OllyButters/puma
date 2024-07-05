#########
# Get the necessary data from Zotero. We can get all the items in the library,
# or just the new ones. The new ones are the ones that are not in the cache.
# It sort of depends how much the data in the Zotero library changes, getting
# the whole lot doesn't seem to take too long, so may be the best option.
######

import logging
from pyzotero import zotero

from config import config
from . import papersCache as pc

################################################################################
# Use the zotero API to get ALL the items in the library.
# This is a bit of a brute force approach, but it does pick up changes that 
# have happened to the zotero library. Querying by updated date might be better.
################################################################################
def getZoteroAll():

    print("Getting all zotero items")
    logging.info('Getting all zotero items')

    zot = zotero.Zotero(config.zotero_id, config.zotero_type, config.zotero_api_key)
    items = zot.everything(zot.top())

    for item in items:
        print('Item: %s | Key: %s' % (item['data']['itemType'], item['data']['key']))
        logging.debug('Item: %s | Key: %s' % (item['data']['itemType'], item['data']['key']))

        # cache as we go
        pc.dumpJson(item['data']['key'], item, 'raw/zotero')


def getZoteroNew():
    print("Getting new zotero items")
    logging.info('Getting new zotero items')

    # Get list of zotero files already in the cache
    zot_cache = pc.getCacheList(filetype='/raw/zotero')

    print("Zotero cache contents:")
    print(zot_cache)
    logging.debug('Zotero cache contents: %s', zot_cache)

    # Get list of zotero items in the online library
    zot = zotero.Zotero(config.zotero_id, config.zotero_type, config.zotero_api_key)
    zot_items = zot.top(format='keys')
    zot_items = zot_items.decode('ascii').split('\n')
    zot_items = list(filter(None, zot_items))

    print("Zotero items:")
    print(zot_items)
    logging.debug('Zotero items: %s', zot_items)

    # Find the missing items
    missing_items = list(set(zot_items) - set(zot_cache))

    print("Missing items:")
    print(missing_items)
    logging.debug('Missing items: %s', missing_items)

    # Get the missing items
    for this_item in missing_items:

        zot_item = zot.top(itemKey = this_item)

        if len(zot_item) == 0:
            print('Item not found: %s' % this_item)
            logging.error('Item not found: %s', this_item)
            continue

        if len(zot_item) > 1:
            print('More than one item found for key: %s' % this_item)
            logging.error('More than one item found for key: %s', this_item)
            continue

        print('Item: %s | Key: %s' % (zot_item[0]['data']['itemType'], zot_item[0]['data']['key']))
        logging.debug('Item: %s | Key: %s' % (zot_item[0]['data']['itemType'], zot_item[0]['data']['key']))

        # cache as we go
        pc.dumpJson(zot_item[0]['data']['key'], zot_item[0], 'raw/zotero')
