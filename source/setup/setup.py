#!/usr/bin/env python2

import os.path
import shutil
import logging
import datetime
import config.config as config


################################################################################
def clean_old_scopus_cache_file():
    # If the force update flag is set in the config then all the scopus data
    # will have been deleted already in the setup phase.

    # Check the age of the exsiting files - scopus doesnt allow old citations
    # so dump the whole file if it is too old.
    cached_scopus_files_list = os.listdir(config.cache_dir + '/raw/scopus/')
    for this_file in cached_scopus_files_list:
        if abs(datetime.datetime.now() - datetime.datetime.fromtimestamp(os.stat(config.cache_dir + '/raw/scopus/' + this_file).st_mtime)) > datetime.timedelta(days=config.scopus_citation_max_age_days):
            print('Deleting: ' + this_file)
            logging.info('Deleting: ' + this_file)
            os.remove(this_file)
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

    html_directories = {"/mesh", "/css", "/papers", "/all_keywords", "/major_keywords", "/country", "/institute", "/metrics", "/wordcloud", "/abstractwordcloud", "/authornetwork", "/help", "/search", "/status"}

    for direct in html_directories:
        if not os.path.exists(config.html_dir + direct):
            os.mkdir(config.html_dir + direct)
