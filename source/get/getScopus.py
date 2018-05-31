#! /usr/bin/env python2

import urllib2
import json
import datetime
import logging
import os

import config.config as config
import papersCache as pc


def cleanScopus(citation_max_life):
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


# Use the elsevier API to get extra info about the paper (including citations).
# Try using the PMID first, if nothing returned then try using the DOI.
def getScopus(zotero_ID, PMID, DOI):

    url = 'http://api.elsevier.com/content/search/scopus'

    try:
        # query scopus with a pmid
        if PMID != "":
            request_string = url + '?apiKey=' + config.scopus_api_key + '&query=PMID(' + str(PMID) + ')'
            logging.info(request_string)
            response = urllib2.urlopen(request_string).read()
            scopus_object = json.loads(response)

            # Check to see if this is just an errorlog
            try:
                error = scopus_object['search-results']['entry'][0]['error']
                if error == 'Result set was empty':
                    logging.info('Result set empty for ' + str(PMID))
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
            if DOI != "":
                try:
                    request_string = url + '?apiKey=' + config.scopus_api_key + '&query=DOI(' + DOI + ')'
                    logging.info(request_string)
                    response = urllib2.urlopen(request_string).read()
                    scopus_object = json.loads(response)

                    # Check to see if this is just an errorlog
                    try:
                        error = scopus_object['search-results']['entry'][0]['error']
                        if error == 'Result set was empty':
                            logging.info('Result set empty for ' + str(DOI))
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
                filename = zotero_ID + '.scopus'
                pc.dumpJson(filename, scopus_object, filetype='raw/scopus')
                return scopus_object
            except:
                print('No scopus ID found.')
                return None
    except:
        return None
