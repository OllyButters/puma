from . import papersZotero as pz
from . import papersCache as pc
from . import getDoi as pd
from . import getPubmed as pm
from . import getScopus as ps
import hashlib
import re
import config.config as config
import logging


# Needs a tidy up and more logging.

def collate():
    # first, check for new papers from zotero repo
    zot = pz.zotPaper()

    # get lists from cache
    zot_cache = pc.getCacheList(filetype='/raw/zotero')
    doi_cache = pc.getCacheList(filetype='/raw/doi')
    pm_cache = pc.getCacheList(filetype='/raw/pubmed')
    scopus_cache = pc.getCacheData(filetype='/raw/scopus')

    zot.collection = config.zotero_collection

    # get list of ALL keys in this zotero library instance
    zot.getPapersKeys()

    new_keys = []
    zotero_papers = []

    # we may want to re-download the data from zotero
    # if config has the 'zotero_get_all' flag set to True, make sure we get all papers not just new ones
    if config.zotero_get_all is True:
        new_keys = zot.papers_keys
    else:
        for num, paper_key in enumerate(zot.papers_keys):
            if paper_key not in zot_cache:
                new_keys.append(paper_key)
            else:
                # get the previously downloaded papers from the cache
                cached_zotero_data = pc.getCacheData(filetype='/raw/zotero', filenames=[paper_key])[paper_key]
                # check itemType - if it's 'note', or 'attachment' we can ignore
                if cached_zotero_data['data']['itemType'] not in ('attachment', 'note'):
                    zotero_papers.append(cached_zotero_data['data'])

    # get all new papers. Note this will not write any to disk until it has got all of them.
    zot.getPapersList(key_list=new_keys)

    # cache the new papers
    for num, paper in enumerate(zot.papers):
        # cache zotero data
        pc.dumpJson(paper['key'], paper, 'raw/zotero')
        # add to new_papers for later doi/pubmed data retrieval
        # check itemType - if it's 'note', or 'attachment' we can ignore
        if paper['data']['itemType'] not in ('attachment', 'note'):
            zotero_papers.append(paper['data'])

    # zotero_papers will have ALL zotero data in it based on cache and newly downloaded data.
    # now check zotero_papers for doi, pubmed id, scopus data and retrieve if required.
    counter = 0
    total_number_zotero_papers = str(len(zotero_papers))

    for paper in zotero_papers:

        # Keep track of how far through we are if watching terminal# check itemType - if it's 'note', or 'attachment' we can ignore
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

        logging.info('\nWorking on '+paper['title']+' (zotero key: '+paper['key']+')')
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
            # doi_filename = hashlib.md5(paper['DOI']).hexdigest()

            logging.info('DOI: ' + doi)
            logging.info('DOI_filename:' + doi_filename)

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
            pmid_matches = re.search(r'PMID: ([0-9]{1,8})', paper['extra'])
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

        # Do scopus data now
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
        logging.info('Merged filename: ' + merged_filename)
        pc.dumpJson(merged_filename, this_merged_paper, 'processed/merged')
