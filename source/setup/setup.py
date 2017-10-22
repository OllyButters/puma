#! /usr/bin/env python

import os.path
import shutil
import config.config as config


###########################################################
# Some of the config settings want to trample cached data,
# this needs to be done elegantly otherwise old data is
# left behind and can get picked up unnecessarily.
def tidy_existing_file_tree(root_dir):

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

    # Merged files. Orphans get left here if deleted from e.g. zotero or
    # the raw cache folder.
    if config.merge_all is True and config.use_cached_merge_only is False:
        if os.path.exists(config.cache_dir + '/processed/merged'):
            shutil.rmtree(config.cache_dir + '/processed/merged')


###########################################################
# Make sure the directory structure is set up first.
# Everything in the cache is grabbed from elsewhere, or built on the fly,
# so it should all be considerd ready to be deleted at any point!
def build_file_tree(root_dir):
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
    if not os.path.exists(root_dir + '/data'):
        os.mkdir(root_dir + '/data')

    if not os.path.exists(config.data_dir):
        os.mkdir(config.data_dir)

    # = Log directory =
    if not os.path.exists(root_dir + '/logs'):
        os.mkdir(root_dir + '/logs')

    # = Output html directory =
    # Html dir
    if not os.path.exists(root_dir + "/html"):
        os.mkdir(root_dir + "/html")

    # Study html dir. This one is inside the previous one.
    if not os.path.exists(config.html_dir):
        os.mkdir(config.html_dir)

    html_directories = {"/mesh", "/css", "/papers", "/all_keywords", "/major_keywords", "/map", "/country", "/city", "/metrics", "/wordcloud", "/abstractwordcloud", "/authornetwork", "/errorlog", "/help", "/search", "status"}

    for direct in html_directories:
        if not os.path.exists(config.html_dir + direct):
            os.mkdir(config.html_dir + direct)
