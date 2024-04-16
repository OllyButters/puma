#!/usr/bin/env python3

import os.path
import shutil
import logging
import datetime
from config import config

# These should be generic functions

################################################################################
def clean_old_zotero_cache_file():
    # Check the age of the existing files - zotero files get updated from time to time
    if os.path.exists(config.cache_dir + '/raw/zotero'):
        logging.debug('Cleaning old zotero files')
        cached_zotero_files_list = os.listdir(config.cache_dir + '/raw/zotero/')
        logging.debug('Found %s  zotero files', str(len(cached_zotero_files_list)))
        for this_file in cached_zotero_files_list:
            if abs(datetime.datetime.now() - datetime.datetime.fromtimestamp(os.stat(config.cache_dir + '/raw/zotero/' + this_file).st_mtime)) > datetime.timedelta(days=config.zotero_data_max_age_days):
                file_path = config.cache_dir + '/raw/zotero/' + this_file
                print('Deleting: ' + file_path)
                logging.info('Deleting: %s', file_path)
                os.remove(file_path)
################################################################################

################################################################################
def clean_old_pubmed_cache_file():
    # Check the age of the existing files - pubmed files get updated from time to time
    if os.path.exists(config.cache_dir + '/raw/pubmed'):
        logging.debug('Cleaning old pubmed files')
        cached_pubmed_files_list = os.listdir(config.cache_dir + '/raw/pubmed/')
        logging.debug('Found %s pubmed files', str(len(cached_pubmed_files_list)))
        for this_file in cached_pubmed_files_list:
            if os.path.isfile(config.cache_dir + '/raw/pubmed/' + this_file): 
                if abs(datetime.datetime.now() - datetime.datetime.fromtimestamp(os.stat(config.cache_dir + '/raw/pubmed/' + this_file).st_mtime)) > datetime.timedelta(days=config.pubmed_data_max_age_days):
                    file_path = config.cache_dir + '/raw/pubmed/' + this_file
                    print('Deleting: ' + file_path)
                    logging.info('Deleting: %s', file_path)
                    os.remove(file_path)
                    file_path = config.cache_dir + '/raw/pubmed/xml/' + this_file + '.xml'
                    os.remove(file_path)
################################################################################

################################################################################
def clean_old_scopus_cache_file():
    # If the force update flag is set in the config then all the scopus data
    # will have been deleted already in the setup phase.

    # Check the age of the existing files - scopus doesnt allow old citations
    # so dump the whole file if it is too old.
    if os.path.exists(config.cache_dir + '/raw/scopus'):
        logging.debug('Cleaning old scopus files')
        cached_scopus_files_list = os.listdir(config.cache_dir + '/raw/scopus/')
        logging.debug('Found %s scopus files', str(len(cached_scopus_files_list)))
        for this_file in cached_scopus_files_list:
            if abs(datetime.datetime.now() - datetime.datetime.fromtimestamp(os.stat(config.cache_dir + '/raw/scopus/' + this_file).st_mtime)) > datetime.timedelta(days=config.scopus_citation_max_age_days):
                file_path = config.cache_dir + '/raw/scopus/' + this_file
                print('Deleting: ' + file_path)
                logging.info('Deleting: %s', file_path)
                os.remove(file_path)
################################################################################


################################################################################
# Some of the config settings want to trample cached data,
# this needs to be done elegantly otherwise old data is
# left behind and can get picked up unnecessarily.
################################################################################
def tidy_existing_file_tree():

    # Raw zotero files. Orphans get left here if deleted from e.g. zotero
    if config.zotero_get_all is True:
        if os.path.exists(config.cache_dir + '/raw/zotero'):
            shutil.rmtree(config.cache_dir + '/raw/zotero')

    # Raw pubmed and DOI files. Orphans get left here if deleted from e.g. zotero
    if config.use_doi_pubmed_cache is False:
        if os.path.exists(config.cache_dir + '/raw/doi'):
            shutil.rmtree(config.cache_dir + '/raw/doi')

        if os.path.exists(config.cache_dir + '/raw/pubmed'):
            shutil.rmtree(config.cache_dir + '/raw/pubmed')

    if config.scopus_force_citation_update is True:
        if os.path.exists(config.cache_dir + '/raw/scopus'):
            shutil.rmtree(config.cache_dir + '/raw/scopus')

    # Merged files. Orphans get left here if deleted from e.g. zotero or
    # the raw cache folder. These get rebuilt everytime now, so trash everything
    # to avoid this.
    # if config.merge_all is True and config.use_cached_merge_only is False:
    if os.path.exists(config.cache_dir + '/processed/merged'):
        shutil.rmtree(config.cache_dir + '/processed/merged')

    if os.path.exists(config.cache_dir + '/processed/cleaned'):
        shutil.rmtree(config.cache_dir + '/processed/cleaned')

    # Delete coverage report if it exists
    if os.path.exists(config.cache_dir + '/coverage_report.html'):
        os.remove(config.cache_dir + '/coverage_report.html')

    if os.path.exists(config.html_dir):
        shutil.rmtree(config.html_dir)


################################################################################
# Make sure the directory structure is set up first.
# Everything in the cache is grabbed from elsewhere, or built on the fly,
# so it should all be considerd ready to be deleted at any point!
################################################################################
def build_file_tree():
    # = Cache directory =
    if not os.path.exists(config.cache_dir):
        os.makedirs(config.cache_dir)

    if not os.path.exists(config.cache_dir + '/raw'):
        os.mkdir(config.cache_dir + '/raw')

    if not os.path.exists(config.cache_dir + '/raw/doi'):
        os.mkdir(config.cache_dir + '/raw/doi')

    if not os.path.exists(config.cache_dir + '/raw/pubmed'):
        os.mkdir(config.cache_dir + '/raw/pubmed')

    if not os.path.exists(config.cache_dir + '/raw/pubmed/xml'):
        os.mkdir(config.cache_dir + '/raw/pubmed/xml')

    if not os.path.exists(config.cache_dir + '/raw/scopus'):
        os.mkdir(config.cache_dir + '/raw/scopus')

    if not os.path.exists(config.cache_dir + '/raw/zotero'):
        os.mkdir(config.cache_dir + '/raw/zotero')

    if not os.path.exists(config.cache_dir + '/processed'):
        os.mkdir(config.cache_dir + '/processed')

    if not os.path.exists(config.cache_dir + '/processed/merged'):
        os.mkdir(config.cache_dir + '/processed/merged')

    if not os.path.exists(config.cache_dir + '/processed/cleaned'):
        os.mkdir(config.cache_dir + '/processed/cleaned')

    if not os.path.exists(config.cache_dir + '/geodata'):
        os.mkdir(config.cache_dir + '/geodata')

    # = Data directory =
    if not os.path.exists(config.data_dir):
        os.makedirs(config.data_dir)

    # = Log directory =
    if not os.path.exists(config.log_dir):
        os.makedirs(config.log_dir)

    # = Output html directory =
    # Html dir
    if not os.path.exists(config.html_dir):
        os.makedirs(config.html_dir)

    html_directories = {"/mesh", "/css", "/papers", "/tags", "/keywords", "/country", "/institute", "/metrics", "/keyword_wordcloud", "/abstractwordcloud", "/authornetwork", "/help", "/search", "/status"}

    for direct in html_directories:
        if not os.path.exists(config.html_dir + direct):
            os.mkdir(config.html_dir + direct)
