import configparser
import sys
import argparse
import os.path

path_to_config_py = os.path.abspath(sys.argv[0])
root_dir = os.path.dirname(os.path.dirname(path_to_config_py))
print('Root directory = ' + root_dir)


# Parse all the config in the settings.ini file and put them into a global variable
# these will be accessible via e.g. config.scopus_api_key in all the other modules
# where they have import config.config as config

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
                        'header_institution_logo_filename': config.get('project_details', 'header_institution_logo_filename')}

    # Scopus settings
    scopus_force_citation_update = config.getboolean('scopus', 'scopus_force_citation_update', fallback = True)
    scopus_citation_max_age_days = int(config.get('scopus', 'scopus_citation_max_age_days', fallback = 7))
    scopus_api_key = config.get('scopus', 'scopus_api_key')

    # Zotero settings
    zotero_id = config.get('zotero_api', 'zotero_id')
    zotero_type = config.get('zotero_api', 'zotero_type')
    zotero_api_key = config.get('zotero_api', 'zotero_api_key')
    zotero_collection = config.get('zotero_api', 'zotero_collection')
    zotero_data_max_age_days = int(config.get('zotero_api', 'zotero_data_max_age_days', fallback = 30))
    zotero_get_all = config.getboolean('zotero_api', 'zotero_get_all', fallback = True)

    # collate settings
    use_doi_pubmed_cache = config.getboolean('collate', 'use_doi_pubmed_cache', fallback = True)

    # Pubmed
    pubmed_email = config.get('pubmed_api', 'pubmed_email')
    pubmed_data_max_age_days = int(config.get('pubmed_api', 'pubmed_data_max_age_days', fallback = 30))

    # logging
    logging_loglevel = config.get('logging', 'loglevel', fallback = 'DEBUG')

    # Pages
    WEB_PAGE_SHOW_INSTITUTE_UK_MAP = config.getboolean('pages', 'web_page_show_institute_UK_map', fallback = True)
    WEB_PAGE_SHOW_INSTITUTE_COUNTRY_MAP = config.getboolean('pages', 'web_page_show_institute_country_map', fallback = True)
    WEB_PAGE_SHOW_ZOTERO_TAGS = config.getboolean('pages', 'web_page_show_zotero_tags', fallback = True)
    web_page_is_in_iframe = config.getboolean('pages', 'web_page_is_in_iframe', fallback = False)
    WEB_PAGE_REPORTS = config.get('pages', 'web_page_reports', fallback = None)
    WEB_PAGE_PROJECT_ABOUT_HTML_FILE = config.get('pages', 'web_page_project_about_html_file', fallback = None)

except Exception as e:
    print('Problem with the settings file')
    print(str(e))
    exit(0)

# Define the file tree names
cache_dir = root_dir + '/cache/' + project_details['short_name']
config_dir = root_dir + '/config'
data_dir = root_dir + '/data/' + project_details['short_name']
html_dir = root_dir + '/html/' + project_details['short_name']
template_dir = root_dir + '/source/web_pages/template'
log_dir = root_dir + '/logs'
