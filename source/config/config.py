#! /usr/bin/env python

import ConfigParser


# Parse all the config in the settings.ini file and put them into a global variable
# these will be accessible via config.scopus_api_key in all the other modules
# where they have import config.config as config
def parse_config(root_dir):
    global scopus_force_citation_update
    global scopus_citation_max_age_days
    global scopus_run_citation
    global scopus_api_key
    global logging_loglevel

    config = ConfigParser.ConfigParser()
    config.read(root_dir + "/config/config.ini_sample")
    print config.sections()
    try:
        # Scopus settings
        scopus_force_citation_update = config.get('scopus', 'scopus_force_citation_update')
        scopus_citation_max_age_days = int(config.get('scopus', 'scopus_citation_max_age_days'))
        scopus_run_citation = config.get('scopus', 'scopus_run_citation')
        scopus_api_key = config.get('scopus', 'scopus_api_key')

        # logging
        logging_loglevel = config.get('logging', 'loglevel')
    except:
        print 'Problem with the settings file'
        exit(0)
