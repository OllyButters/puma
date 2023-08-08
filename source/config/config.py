#!/usr/bin/env python3
import configparser
import sys
import argparse
import os.path


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

    global metrics_study_start_year
    global metrics_study_current_year

    global zotero_id
    global zotero_type
    global zotero_api_key
    global zotero_collection

    # collation settings
    global zotero_get_all
    global use_doi_pubmed_cache

    global pubmed_email

    # Some data cannot be dispayed on a public facing website - e.g. all the data in merged json file.
    global web_page_public_facing
    global web_page_is_in_iframe
    global web_page_show_institute_country_map
    global web_page_show_author_network
    global web_page_show_zotero_tags

    global network_create_networks

    # See if there is a config file specified from the command line, if not then
    # use the default config file name.
    config_file_name = 'config.ini'

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="Specify a config file. Just want the name not the path.")
    args = parser.parse_args()
    if args.config:
        print('Config file set to: '+args.config)
        config_file_name = args.config

    # Check config file exists
    config_file_path = root_dir + "/config/" + config_file_name
    if os.path.isfile(config_file_path):
        print('Config file ok!: ' + config_file_path)
    else:
        print('Config file doesnt exist!: ' + config_file_path)
        sys.exit()

    config = configparser.ConfigParser()
    config.read(config_file_path)
    try:
        # Project Details
        project_details = {'name': config.get('project_details', 'name'),
                           'short_name': config.get('project_details', 'short_name'),
                           'colour_hex_primary': config.get('project_details', 'colour_hex_primary'),
                           'colour_hex_secondary': config.get('project_details', 'colour_hex_secondary'),
                           'header_institution_url': config.get('project_details', 'header_institution_url'),
                           'header_institution_name': config.get('project_details', 'header_institution_name'),
                           'header_institution_logo_filename': config.get('project_details', 'header_institution_logo_filename'),
                           'side_image_filename': config.get('project_details', 'side_image_filename'),
                           'side_image_link': config.get('project_details', 'side_image_link')}

        # Scopus settings
        scopus_force_citation_update = config.get('scopus', 'scopus_force_citation_update')
        scopus_citation_max_age_days = int(config.get('scopus', 'scopus_citation_max_age_days'))
        scopus_run_citation = config.get('scopus', 'scopus_run_citation')
        scopus_api_key = config.get('scopus', 'scopus_api_key')

        # Zotero settings
        zotero_id = config.get('zotero_api', 'zotero_id')
        zotero_type = config.get('zotero_api', 'zotero_type')
        zotero_api_key = config.get('zotero_api', 'zotero_api_key')
        zotero_collection = config.get('zotero_api', 'zotero_collection')

        # collate settings
        zotero_get_all = config.getboolean('collate', 'zotero_get_all')
        use_doi_pubmed_cache = config.getboolean('collate', 'use_doi_pubmed_cache')

        # Pubmed
        pubmed_email = config.get('pubmed_api', 'pubmed_email')

        # logging
        logging_loglevel = config.get('logging', 'loglevel')

        # Metrics
        metrics_study_start_year = int(config.get('metrics', 'metrics_study_start_year'))
        metrics_study_current_year = int(config.get('metrics', 'metrics_study_current_year'))

        # Pages
        web_page_show_author_network = config.getboolean('pages', 'web_page_show_author_network')
        web_page_show_institute_country_map = config.getboolean('pages', 'web_page_show_institute_country_map')
        web_page_show_zotero_tags = config.getboolean('pages', 'web_page_show_zotero_tags')
        web_page_is_in_iframe = config.getboolean('pages', 'web_page_is_in_iframe')

        # Networks
        network_create_networks = config.getboolean('networks', 'create_networks')

    except Exception as e:
        print('Problem with the settings file')
        print(str(e))
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
    data_dir = root_dir + '/data/' + project_details['short_name']
    html_dir = root_dir + '/html/' + project_details['short_name']
    template_dir = root_dir + '/source/web_pages/template'
    log_dir = root_dir + '/logs'
