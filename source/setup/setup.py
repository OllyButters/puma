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

    if not os.path.exists(config.cache_dir):
        os.mkdir(config.cache_dir)

    if not os.path.exists(config.cache_dir + '/geodata'):
        os.mkdir(config.cache_dir + '/geodata')

    if not os.path.exists(root_dir + '/data'):
        os.mkdir(root_dir + '/data')

    # Log directory
    if not os.path.exists(root_dir + '/logs'):
        os.mkdir(root_dir + '/logs')

    # Output html
    if not os.path.exists(root_dir + "/html"):
        os.mkdir(root_dir + "/html")

    if not os.path.exists(config.html_dir):
        os.mkdir(config.html_dir)

    if not os.path.exists(config.html_dir + '/mesh'):
        os.mkdir(config.html_dir + '/mesh')

    if not os.path.exists(config.html_dir + '/css'):
        os.mkdir(config.html_dir + '/css')

    if not os.path.exists(config.html_dir + '/papers'):
        os.mkdir(config.html_dir + '/papers')

    if not os.path.exists(config.html_dir + '/all_keywords'):
        os.mkdir(config.html_dir + '/all_keywords')

    if not os.path.exists(config.html_dir + '/major_keywords'):
        os.mkdir(config.html_dir + '/major_keywords')

    if not os.path.exists(config.html_dir + '/map'):
        os.mkdir(config.html_dir + '/map')

    if not os.path.exists(config.html_dir + '/country'):
        os.mkdir(config.html_dir + '/country')

    if not os.path.exists(config.html_dir + '/city'):
        os.mkdir(config.html_dir + '/city')

    if not os.path.exists(config.html_dir + '/metrics'):
        os.mkdir(config.html_dir + '/metrics')

    if not os.path.exists(config.html_dir + '/wordcloud'):
        os.mkdir(config.html_dir + '/wordcloud')

    if not os.path.exists(config.html_dir + '/abstractwordcloud'):
        os.mkdir(config.html_dir + '/abstractwordcloud')

    if not os.path.exists(config.html_dir + '/authornetwork'):
        os.mkdir(config.html_dir + '/authornetwork')

    if not os.path.exists(config.html_dir + '/errorlog'):
        os.mkdir(config.html_dir + '/errorlog')
