import papersZotero as pz
import papersCache as pc
import getDoi as pd
import getPubmed as pm
import papersMerge as pMerge
import hashlib
import re
import copy
import json
import config.config as config
import logging


def collate():
    # first, check if config.use_cached_merged_only is set to True
    # if so, just load these and return
    if config.use_cached_merge_only is True:
        # get list of currently merged papers
        merged_list = pc.getCacheList(filetype='/processed/merged')
        merged_papers = []
        logging.info("COLLATE: use_cached_merge_only set to 1 so loading papers straight from processed/merged")
        print "use_cached_merge_only set to 1 so loading papers straight from processed/merged"
        for merged_paper in merged_list:
            logging.info("Loading cached merged paper: "+merged_paper)
            print "Loading cached merged paper: "+merged_paper+"."
            merged_papers.append(pc.getCacheData(filetype='/processed/merged', filenames=[merged_paper])[merged_paper])

        return merged_papers

    # first, check for new papers from zotero repo
    zot = pz.zotPaper()

    # get list from cache
    zot_cache = pc.getCacheList(filetype='/raw/zotero')
    doi_cache = pc.getCacheList(filetype='/raw/doi')
    pm_cache = pc.getCacheList(filetype='/raw/pubmed')

    zot.collection = config.zotero_collection

    # get list of all keys in this zotero instance
    zot.getPapersKeys()

    new_keys = []
    new_papers = []

    # we may want to re-download the data from zotero
    # if config has the 'zotero_get_all' flag set to True, make sure we get all papers not just new ones (i.e. load from cached file)
    if config.zotero_get_all is True:
        new_keys = zot.papers_keys
    else:
        for num, paper_key in enumerate(zot.papers_keys):
            if paper_key not in zot_cache:
                new_keys.append(paper_key)
            else:
                # get the previously downloaded papers from the cache
                new_paper = pc.getCacheData(filetype='/raw/zotero', filenames=[paper_key])[paper_key]
                # check itemType - if it's 'note', we can ignore
                # if new_paper['data']['itemType'] != 'note':
                if new_paper['data']['itemType'] not in ('attachment', 'note'):
                    new_papers.append(new_paper['data'])

    # get all new papers
    zot.getPapersList(key_list=new_keys)

    # cache the new papers
    for num, paper in enumerate(zot.papers):
        # cache zotero data
        pc.dumpJson(paper['key'], paper, 'raw/zotero')
        # add to new_papers for later doi/pubmed data retrieval
        # check itemType - if it's 'note', we can ignore
        if paper['data']['itemType'] != 'note':
            new_papers.append(paper['data'])

    # now check new_papers for doi or pubmed id and retrieve if required
    for paper in new_papers:
        # add an IDs key to paper and populate with zotero key
        # plus empty keys for doi, pmid and hash
        paper['IDs'] = {
          'zotero': paper['key'],
          'DOI': '',
          'DOI_filename': '',
          'PMID': '',
          'hash': '',
        }

        logging.info('\nWorking on '+paper['title']+' (zotero key: '+paper['key']+')')
        logging.info('Getting doi/pubmed data.')
        print 'Getting doi/pubmed data for paper: '+paper['title']+' (zotero key: '+paper['key']+')'

        paper['doi_data'] = {}
        paper['pmid_data'] = {}
        if 'DOI' in paper and paper['DOI'] != "":
            # check if this is the full url, and if so, strip the parent
            # note that this is done again in getDoi.getDoi, but will not find the filename in the cache if we don't process this
            check_doi = re.match(r'^http://dx\.doi\.org/', paper['DOI'])
            if check_doi is not None:
                doi = re.sub(r'^http://dx\.doi\.org/', '', doi)

            # as doi's use '/' chars, we do an md5 of the doi as the filename
            print ((paper['DOI']).encode("ascii", "ignore"))
            doi_filename = hashlib.md5((paper['DOI']).encode("ascii", "ignore")).hexdigest()
            # doi_filename = hashlib.md5(paper['DOI']).hexdigest()

            logging.info('DOI: ' + paper['DOI'])
            logging.info('DOI_filename:' + doi_filename)

            # add these ids to the IDs dict
            paper['IDs']['DOI'] = paper['DOI']
            paper['IDs']['DOI_filename'] = doi_filename

            # check if paper data in doi cache (only if config.use_pubmed_doi_cache is not True
            if doi_filename not in doi_cache or config.use_doi_pubmed_cache is not True:
                logging.debug('Downloading (or redownloading) DOI data.')
                doi_paper = pd.getDoi(paper['DOI'])
                paper['doi_data'] = doi_paper
                # data is automatically cached by getDoi
            else:
                logging.debug('Getting cached DOI data.')
                paper['doi_data'] = pc.getCacheData(filetype='/raw/doi', filenames=[doi_filename])[doi_filename]

        # get pubmed data
        if 'extra' in paper:
            pmid_matches = re.search(r'PMID: ([0-9]{1,8})', paper['extra'])
            if pmid_matches is not None:
                paper['pmid'] = pmid_matches.group(1)
                # add the pmid to the IDs dict
                paper['IDs']['PMID'] = paper['pmid']

                # check if paper data in pm cache (only if config.use_doi_pubmed_cache is not True)
                if paper['pmid'] not in pm_cache or config.use_doi_pubmed_cache is not True:
                    pm_paper = pm.getPubmed(paper['pmid'])
                    paper['pmid_data'] = pm_paper
                    # data is automatically cached by getPubmed
                else:
                    paper['pmid_data'] = pc.getCacheData(filetype='/raw/pubmed', filenames=[paper['pmid']])[paper['pmid']]

    # now do merge data
    merged_papers = {}

    # get list of currently merged papers
    merged_list = pc.getCacheList(filetype='/processed/merged')
    # merged_papers = []
    template_file = open(config.config_dir+'/data-doi-template', 'r')
    template = json.load(template_file)
    template_file.close()
    # template = pc.getCacheData(filenames=['data-doi-template'])['data-doi-template']

    logging.info('Starting merging process.')

    for paper in new_papers:
        logging.debug('Merging ' + paper['IDs']['zotero'])
        merged_paper = {}
        doi_data = paper['doi_data']
        pmid_data = paper['pmid_data']
        # delete these from paper
        del paper['doi_data']
        del paper['pmid_data']

        # temporarily remove the IDs dict (for later replacement)
        ids = paper['IDs']
        del paper['IDs']

        # check if title is empty
        # if so, check doi_data and pmid_data for title (in that order)
        # and use instead
        if paper['title'] == '':
            logging.info('title is blank, trying to hack something together.')
            try:
                if doi_data['title'] != '' and doi_data['title'] is not None:
                    paper['title'] = doi_data['title']
            except (KeyError, TypeError):
                try:
                    if pmid_data['MedlineCitation']['Article']['ArticleTitle'] != '' and pmid_data['MedlineCitation']['Article']['ArticleTitle'] is not None:
                        paper['title'] = pmid_data['MedlineCitation']['Article']['ArticleTitle']
                except (KeyError, TypeError):
                    logging.warn('No title for Zotero id: ' + paper['key'])
                    pass

        # create new filename (md5 of title)
        filename = hashlib.md5(paper['title'].encode('ascii', 'ignore')).hexdigest()
        # make sure we add the filename to ids
        ids['hash'] = filename

        # if we aren't set to merge all papers, ignore existing files
        if config.merge_all is False:
            if filename in merged_list:
                logging.info("Merged file: "+filename+" already exists. Ignoring as merge_all not set to 1 in config.ini. File being loaded from cache.")
                print "Merged file: "+filename+" already exists. Ignoring as merge_all not set to 1 in config.ini. File being loaded from cache."
                paper = pc.getCacheData(filetype='/processed/merged', filenames=[filename])[filename]
                continue

        logging.info('Merged filename: ' + filename)
        print "Merging to filename: "+filename

        # the merge process uses the merge class and jsonpath mappings (see config dir)
        # data is merged separately from each source to the blank template (data-doi-template in config)
        # these separate 'doi-formatted' datasets are then merged together (pubmed into zotero, doi into pubmed/zotero)

        mgr = pMerge.Merge()

        ###############
        # Pubmed first. This goes into 'pmid_data' section of merged_paper
        map_file = open(config.config_dir+'/data-pubmed-doi-jsonpath', 'r')
        mgr.mapping = json.load(map_file)
        map_file.close()
        # mgr.mapping = pc.getCacheData(filenames=['data-pubmed-doi-jsonpath'])['data-pubmed-doi-jsonpath']
        mgr.src = pmid_data
        mgr.dest = copy.deepcopy(template)
        logging.info('Copying pubmed data.')
        print 'Merging pubmed'
        mgr.mapSrc()
        merged_paper['pmid_data'] = copy.deepcopy(mgr.dest)

        #################
        # now do zotero data. This goes into 'zotero_data' section of merged_paper
        map_file = open(config.config_dir+'/data-zotero-doi-jsonpath', 'r')
        mgr.mapping = json.load(map_file)
        map_file.close()
        # mgr.mapping = pc.getCacheData(filenames=['data-zotero-doi-jsonpath'])['data-zotero-doi-jsonpath']
        mgr.src = paper
        mgr.dest = copy.deepcopy(template)
        logging.info('Copying zotero data.')
        print 'Merging zotero'
        mgr.mapSrc()
        merged_paper['zotero_data'] = copy.deepcopy(mgr.dest)

        #################
        # now merge doi_data, merged_paper
        # starting with zotero_data, we merge using an empty mapping and setting the dest as zotero_data and src as pmid_data
        # This will mean the pmid_data will overwrite zotero_data where there are overlapping fields.
        # If there are gaps in the zotero data then this will hopefully fill them in.
        mgr.mapping = {}
        mgr.src = merged_paper['pmid_data']
        mgr.dest = merged_paper['zotero_data']
        logging.info('Merging (overwriting) pubmed data onto zotero.')
        print 'Merging pubmed into zotero'
        mgr.mapSrc()
        merged_paper['pmid_zotero_data'] = copy.deepcopy(mgr.dest)

        ##################
        # now set the src to be doi_data and merge to the template (otherwise the output data gets oddly formatted)
        mgr.mapping = {}
        mgr.src = doi_data
        mgr.dest = copy.deepcopy(template)
        logging.info('Copying DOI data.')
        print 'Merging doi'
        mgr.mapSrc()
        merged_paper['doi_data'] = copy.deepcopy(mgr.dest)

        #################
        # now set the src to be doi_data and merge to the dest (pmid/zotero data)
        mgr.mapping = {"$.DOI": "$.DOI", "$.title": "$.title"}
        mgr.dest = merged_paper['doi_data']
        mgr.src = merged_paper['pmid_zotero_data']
        logging.info('Merging (overwriting) DOI data onto already merged zotero/pmid data.')
        print 'Merging doi into pubmed/zotero'
        mgr.mapSrc()

        # now do some restructuring do make sure we dump data in the correct
        # object format for later processing -
        # {
        #   'IDs': [dict of all ids],
        #   'merged': [the final merged object (copy of mgr.dest)]
        # }
        merged_paper = {
          'IDs': ids,
          'merged': copy.deepcopy(mgr.dest)
        }
        merged_papers[filename] = merged_paper
        # merged_papers[filename] = copy.deepcopy(mgr.dest)
        pc.dumpJson(filename, merged_papers[filename], 'processed/merged')

    return merged_papers

if __name__ == "__main__":
    print "Collate data"

    collate()
