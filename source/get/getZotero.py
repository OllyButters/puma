import hashlib
import re
import logging
from pyzotero import zotero

from config import config as config
from . import papersCache as pc

#config.build_config_variables("../")

allowed_item_types = ['journalArticle']

################################################################################
# Use the zotero API to get ALL the items in the library.
# This is a bit of a brute force approach, but it does pick up changes that 
# have happened to the zotero library. Querying by updated date might be better.
################################################################################
def getZoteroAll():

    zot = zotero.Zotero(config.zotero_id, config.zotero_type, config.zotero_api_key)
    items = zot.top(limit=50)

    for item in items:
        print('Item: %s | Key: %s' % (item['data']['itemType'], item['data']['key']))
        # cache as we go
        pc.dumpJson(item['data']['key'], item, 'raw/zotero')

