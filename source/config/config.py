#! /usr/bin/env python
import ConfigParser
import os
import sys

# Parse all the config in the settings.ini file and put them into a global variable
# these will be accessible via config.scopus_api_key in all the other modules
# where they have import config.config as config
def build_config_variables(root_dir):

    # Parse the config file
    global project_details

    global scopus_force_citation_update
    global scopus_citation_max_age_days
    global scopus_run_citation
    global scopus_api_key

    global logging_loglevel

    global google_maps_api_key

    global metrics_study_start_year
    global metrics_study_current_year

    global zotero_id
    global zotero_type
    global zotero_api_key
    global zotero_collection

    global pubmed_email

    global page_show_author_network

    config = ConfigParser.ConfigParser()
    config.read(root_dir + "/config/config.ini")
    try:
        # Project Details
        project_details = {'name': config.get('project_details', 'name'), 'short_name': config.get('project_details', 'short_name'), 'colour_hex_primary': config.get('project_details', 'colour_hex_primary'), 'colour_hex_secondary': config.get('project_details', 'colour_hex_secondary'), 'header_image_url': config.get('project_details', 'header_image_url'), 'header_institution': config.get('project_details', 'header_institution'), 'header_institution_url': config.get('project_details', 'header_institution_url'), 'side_image_url': config.get('project_details', 'side_image_url'), 'side_image_link': config.get('project_details', 'side_image_link')}

        # Scopus settings
        scopus_force_citation_update = config.get('scopus', 'scopus_force_citation_update')
        scopus_citation_max_age_days = int(config.get('scopus', 'scopus_citation_max_age_days'))
        scopus_run_citation = config.get('scopus', 'scopus_run_citation')
        scopus_api_key = config.get('scopus', 'scopus_api_key')

        #Zotero settings
        zotero_id = config.get('zotero_api', 'zotero_id')
        zotero_type = config.get('zotero_api', 'zotero_type')
        zotero_api_key = config.get('zotero_api', 'zotero_api_key')
        zotero_collection = config.get('zotero_api', 'zotero_collection')

        #Pubmed
        pubmed_email = config.get('pubmed_api', 'pubmed_email')

        # logging
        logging_loglevel = config.get('logging', 'loglevel')

        # Mapping
        google_maps_api_key = config.get('google_apis', 'google_maps_api_key')

        # Metrics
        metrics_study_start_year = int(config.get('metrics', 'metrics_study_start_year'))
        metrics_study_current_year = int(config.get('metrics', 'metrics_study_current_year'))

        # Mages
        page_show_author_network = config.get('pages', 'page_show_author_network')

    except:
        print 'Problem with the settings file'
        exit(0)

    # Define the file tree names
    global cache_dir
    global config_dir
    global data_dir
    global html_dir
    global template_dir
    global log_dir

    cache_dir = root_dir + '/cache/' + project_details['short_name']
    config_dir = root_dir + '/config'
    data_dir = root_dir + '/data'
    html_dir = root_dir + '/html/' + project_details['short_name']
    template_dir = root_dir + '/source/html/template'
    log_dir = root_dir + '/logs'
