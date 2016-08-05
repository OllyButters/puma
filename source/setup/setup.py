#! /usr/bin/env python

import os.path
import config.config as config


###########################################################
# Make sure the directory structure is set up first.
# Everything in the cache is grabbed from elsewhere, or built on the fly,
# so it should all be considerd ready to be deleted at any point!
def build_file_tree(root_dir):
    # Cache first
    # if (os.path.exists(root_dir + '/cache') is False):
    #    os.mkdir(root_dir + '/cache')

    if (os.path.exists(config.cache_dir) is False):
        os.mkdir(config.cache_dir)

    if (os.path.exists(config.cache_dir + '/geodata') is False):
        os.mkdir(config.cache_dir + '/geodata')

    if (os.path.exists(root_dir + '/data') is False):
        os.mkdir(root_dir + '/data')

    # Log directory
    if (os.path.exists(root_dir + '/logs') is False):
        os.mkdir(root_dir + '/logs')

    # Output html
    if not os.path.exists(root_dir + '/html'):
        os.mkdir(root_dir + '/html')

    if not os.path.exists(root_dir + '/html/mesh'):
        os.mkdir(root_dir + '/html/mesh')

    if not os.path.exists(root_dir + '/html/css'):
        os.mkdir(root_dir + '/html/css')

    if not os.path.exists(root_dir + '/html/papers'):
        os.mkdir(root_dir + '/html/papers')

    if not os.path.exists(root_dir + '/html/all_keywords'):
        os.mkdir(root_dir + '/html/all_keywords')

    if not os.path.exists(root_dir + '/html/major_keywords'):
        os.mkdir(root_dir + '/html/major_keywords')

    if not os.path.exists(root_dir + '/html/map'):
        os.mkdir(root_dir + '/html/map')

    if not os.path.exists(root_dir + '/html/country'):
        os.mkdir(root_dir + '/html/country')

    if not os.path.exists(root_dir + '/html/city'):
        os.mkdir(root_dir + '/html/city')

    if not os.path.exists(root_dir + '/html/metrics'):
        os.mkdir(root_dir + '/html/metrics')

    if not os.path.exists(root_dir + '/html/wordcloud'):
        os.mkdir(root_dir + '/html/wordcloud')

    if not os.path.exists(root_dir + '/html/abstractwordcloud'):
        os.mkdir(root_dir + '/html/abstractwordcloud')

    if not os.path.exists(root_dir + '/html/authornetwork'):
        os.mkdir(root_dir + '/html/authornetwork')

    if not os.path.exists(root_dir + '/html/errorlog'):
        os.mkdir(root_dir + '/html/errorlog')
