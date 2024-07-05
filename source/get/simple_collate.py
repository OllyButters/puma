import hashlib
import re
import logging

from config import config
from . import papersCache as pc
from . import getDoi as pd
from . import getPubmed as pm
from . import getScopus as ps
from . import getZotero as gz

allowed_item_types = ['journalArticle']

# Needs a tidy up and more logging.

def collate():


    # Rebuild zotero cache as required
    if config.zotero_get_all is True:
        gz.getZoteroAll()
    #elif config.zotero_get_updated is True:
        # get items which have been MODIFIED
    #elif config.zotero_get_aged is True:
        # replace old aged items on disk.
    else:
        # get items which are new
        gz.getZoteroNew()

    # get lists from cache
    zot_cache = pc.getCacheList(filetype='/raw/zotero')
    doi_cache = pc.getCacheList(filetype='/raw/doi')
    pm_cache = pc.getCacheList(filetype='/raw/pubmed')
    scopus_cache = pc.getCacheData(filetype='/raw/scopus')

    zotero_papers = []
        
    # Read the data in from the cache
    for this_file in zot_cache:
        cached_zotero_data = pc.getCacheData(filetype='/raw/zotero', filenames=[this_file])[this_file]
        if cached_zotero_data['data']['itemType'] in allowed_item_types:
            zotero_papers.append(cached_zotero_data['data'])

    # zotero_papers will have ALL zotero data in it based on cache and newly downloaded data.
    # now check zotero_papers for doi, pubmed id, scopus data and retrieve if required.
    counter = 0
    total_number_zotero_papers = str(len(zotero_papers))
    scopus_quota_reached = False

    for paper in zotero_papers:

        # Keep track of how far through we are if watching terminal
        counter = counter + 1
        print('\nOn ' + str(counter) + '/' + total_number_zotero_papers)

        this_merged_paper = {}

        # add an IDs key to paper and populate with zotero key
        # plus empty keys for doi, pmid and hash
        this_merged_paper['IDs'] = {
          'zotero': paper['key'],
          'DOI': '',
          'DOI_filename': '',
          'PMID': '',
          'scopus': '',
          'hash': '',
        }

        logging.info('\nWorking on %s (zotero key: %s)', paper['title'], paper['key'])
        logging.info('Getting doi/pubmed/scopus data.')
        print('Getting doi/pubmed/scopus data for paper: ' + paper['title'] + ' (zotero key: ' + paper['key'] + ')')

        # Couple of placeholders
        this_merged_paper['raw'] = {}
        this_merged_paper['raw']['zotero_data'] = paper
        this_merged_paper['raw']['doi_data'] = {}
        this_merged_paper['raw']['pmid_data'] = {}
        this_merged_paper['raw']['scopus_data'] = {}

        if 'DOI' in paper and paper['DOI'] != "":
            # check if this is the full url (it could have been entered into zotero
            # with full url), and if so, strip the parent.
            # note that this is done again in getDoi.getDoi, but will not find the filename in the cache if we don't process this
            check_doi = re.match(r'^https?://(dx\.)?doi\.org/', paper['DOI'])
            if check_doi is not None:
                doi = re.sub(r'^https?://(dx\.)?doi\.org/', '', paper['DOI'])
            else:
                doi = paper['DOI']

            # as dois use '/' chars, we do an md5 of the doi as the filename
            print('DOI:' + doi)
            doi_filename = hashlib.md5((doi).encode("ascii", "ignore")).hexdigest()

            logging.info('DOI: %s', doi)
            logging.info('DOI_filename: %s', doi_filename)

            # add these ids to the IDs dict
            this_merged_paper['IDs']['DOI'] = doi
            this_merged_paper['IDs']['DOI_filename'] = doi_filename

            # check if paper data in doi cache (only if config.use_pubmed_doi_cache is not True
            if doi_filename not in doi_cache or config.use_doi_pubmed_cache is not True:
                print('DOI: Downloading.')
                logging.debug('Downloading (or redownloading) DOI data.')
                doi_paper = pd.getDoi(doi)
                this_merged_paper['raw']['doi_data'] = doi_paper
                print('DOI: Success')
                # data is automatically cached by getDoi
            else:
                print('DOI: Getting from cache.')
                logging.debug('Getting cached DOI data.')
                this_merged_paper['raw']['doi_data'] = pc.getCacheData(filetype='/raw/doi', filenames=[doi_filename])[doi_filename]

        # get pubmed data
        if 'extra' in paper:
            pmid_matches = re.search(r'PMID\s*:\s*([0-9]{1,8})', paper['extra'])
            if pmid_matches is not None:
                paper['pmid'] = pmid_matches.group(1)
                print('PMID: ' + paper['pmid'])
                # add the pmid to the IDs dict
                this_merged_paper['IDs']['PMID'] = paper['pmid']

                # check if paper data in pm cache (only if config.use_doi_pubmed_cache is not True)
                if paper['pmid'] not in pm_cache or config.use_doi_pubmed_cache is not True:
                    print('PMID: Downloading.')
                    logging.debug('Downloading (or redownloading) PMID data.')
                    pm_paper = pm.getPubmed(paper['pmid'])
                    this_merged_paper['raw']['pmid_data'] = pm_paper
                    print('PMID: Success')
                    # data is automatically cached by getPubmed
                else:
                    print('PMID: Getting from cache.')
                    this_merged_paper['raw']['pmid_data'] = pc.getCacheData(filetype='/raw/pubmed', filenames=[paper['pmid']])[paper['pmid']]
            else:
                logging.debug("No PMID found in: %s", paper['extra'])

        # Do scopus data now
        if scopus_quota_reached:
            print("Skipping Scopus as API quota reached")
            logging.warning("Skipping Scopus as API quota reached")
        else:
            if this_merged_paper['IDs']['PMID'] != '' or this_merged_paper['IDs']['DOI'] != '':
                # check if paper data in scopus cache
                scopus_cache_filename = this_merged_paper['IDs']['zotero'] + '.scopus'
                if scopus_cache_filename not in scopus_cache:
                    print('Scopus: Downloading.')
                    logging.debug('Downloading (or redownloading) Scopus data.')
                    logging.debug(this_merged_paper['IDs']['zotero'])
                    logging.debug(this_merged_paper['IDs']['PMID'])
                    logging.debug(this_merged_paper['IDs']['DOI'])
                    scopus_paper = ps.getScopus(this_merged_paper['IDs']['zotero'], this_merged_paper['IDs']['PMID'], this_merged_paper['IDs']['DOI'])

                    # If the scopus quota is hit then stop querying it
                    if scopus_paper == "QUOTAEXCEEDED":
                        scopus_quota_reached = True
                        print("Scopus API quota exceeded. No more Scopus queries for this run.")
                        logging.info("Scopus API quota exceeded. No more Scopus queries for this run.")
                        continue

                    this_merged_paper['raw']['scopus_data'] = scopus_paper
                    print('Scopus: Success')
                    # data is automatically cached by getScopus
                else:
                    print('Scopus: Getting from cache.')
                    this_merged_paper['raw']['scopus_data'] = pc.getCacheData(filetype='/raw/scopus', filenames=scopus_cache_filename)[scopus_cache_filename]

                try:
                    # Lets grab the scopus ID while we are here
                    this_merged_paper['IDs']['scopus'] = this_merged_paper['raw']['scopus_data']['search-results']['entry'][0]['eid']
                except:
                    pass

        # May as well just write out the whole merged file everytime.
        merged_filename = this_merged_paper['IDs']['zotero'] + '.merged'
        logging.info('Merged filename: %s', merged_filename)
        pc.dumpJson(merged_filename, this_merged_paper, 'processed/merged')
