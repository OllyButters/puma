#! /usr/bin/env python2

import urllib2
import json
import datetime
import logging
import os

import config.config as config
import papersCache as pc


# Use the elsevier API to get extra info about the paper (including citations).
# Try using the PMID first, if nothing returned then try using the DOI.
def getScopus(papers, api_key, citation_max_life, force_update, error_log):

    url = 'http://api.elsevier.com/content/search/scopus'

    print('Getting Scopus data')

    # If the force update flag is set in the config then all the scopus data
    # will have been deleted already in the setup phase.

    # Check the age of the exsiting files - scopus doesnt allow old citations
    # so dump the whole file if it is too old.
    print(datetime.datetime.now())
    cached_scopus_files_list = os.listdir(config.cache_dir + '/raw/scopus/')
    for this_file in cached_scopus_files_list:
        print(os.stat(config.cache_dir + '/raw/scopus/' + this_file).st_mtime)
        print(datetime.datetime.fromtimestamp(os.stat(config.cache_dir + '/raw/scopus/' + this_file).st_mtime))
        if abs(datetime.datetime.now() - datetime.datetime.fromtimestamp(os.stat(config.cache_dir + '/raw/scopus/' + this_file).st_mtime)) > datetime.timedelta(days=citation_max_life):
            print('Will delete')

    number_papers_to_process = len(papers)
    counter = 0
    for this_paper in papers:
        counter = counter + 1
        logging.info('\nScopus on # ' + str(counter) + ' of ' + str(number_papers_to_process) + ' ' + this_paper['IDs']['zotero'])
        print('\nScopus on # ' + str(counter) + ' of ' + str(number_papers_to_process) + ' ' + this_paper['IDs']['zotero'])

        # Check if this file already exists, if it does then skip it.
        filename = config.cache_dir + '/raw/scopus/' + this_paper['IDs']['zotero'] + '.scopus'
        if os.path.isfile(filename):
            logging.info(filename + ' alrady exists, skipping.')
            continue

        try:
            # query scopus with a pmid
            if this_paper['IDs']['PMID'] != "":
                request_string = url + '?apiKey=' + api_key + '&query=PMID(' + this_paper['IDs']['PMID'] + ')'
                logging.info(request_string)
                response = urllib2.urlopen(request_string).read()
                scopus_object = json.loads(response)

                # Check to see if this is just an errorlog
                try:
                    error = scopus_object['search-results']['entry'][0]['error']
                    if error == 'Result set was empty':
                        logging.info('Result set empty for ' + str(this_paper['IDs']['PMID']))
                        del(scopus_object)
                except:
                    pass
                if scopus_object:
                    logging.info('Scopus data got via PMID.')
                else:
                    raise

        except:
            # querying with PMID failed, so try DOI.
            try:
                if this_paper['IDs']['DOI'] != "":
                    try:
                        request_string = url + '?apiKey=' + api_key + '&query=DOI(' + this_paper['IDs']['DOI'] + ')'
                        logging.info(request_string)
                        response = urllib2.urlopen(request_string).read()
                        scopus_object = json.loads(response)

                        # Check to see if this is just an errorlog
                        try:
                            error = scopus_object['search-results']['entry'][0]['error']
                            if error == 'Result set was empty':
                                logging.info('Result set empty for ' + str(this_paper['IDs']['PMID']))
                                del(scopus_object)
                        except:
                            pass
                        if scopus_object:
                            logging.info('Scopus data got via DOI.')
                    except:
                        print('Didnt work with DOI either.')
                        logging.info('Scopus data failed to get.')
            except:
                # Nothing found for this papers
                pass

        try:
            # Hopefully have scopus object now.
            if scopus_object != "":
                # Require a scopus ID
                try:
                    # scopus_id = scopus_object['search-results']['entry'][0]['eid']
                    # cannot guarantee format of scopus ID so do an md5 of it for the filename
                    # filename = hashlib.md5(scopus_id).hexdigest()
                    filename = this_paper['IDs']['zotero'] + '.scopus'
                    pc.dumpJson(filename, scopus_object, filetype='raw/scopus')
                except:
                    print('No scopus ID found.')
        except:
            pass

    exit(1)
