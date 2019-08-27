#!/usr/bin/env python2

import urllib2
import json
import logging

import config.config as config
import papersCache as pc


################################################################################
# Use the elsevier API to get extra info about the paper (including citations).
# Try using the PMID first, if nothing returned then try using the DOI.
################################################################################
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
                    logging.warn(scopus_object)
                    del(scopus_object)
            except:
                pass
            if scopus_object:
                logging.info('Scopus data got via PMID.')
            else:
                raise
        else:
            raise
    except:
        # querying with PMID failed, so try DOI.
        try:
            if DOI != "":
                try:
                    logging.debug('now here')
                    request_string = url + '?apiKey=' + config.scopus_api_key + '&query=DOI(' + DOI + ')'
                    logging.info(request_string)
                    response = urllib2.urlopen(request_string).read()
                    scopus_object = json.loads(response)

                    # Check to see if this is just an errorlog
                    try:
                        error = scopus_object['search-results']['entry'][0]['error']
                        if error == 'Result set was empty':
                            logging.info('Result set empty for ' + str(DOI))
                            logging.warn(scopus_object)
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
                filename = zotero_ID + '.scopus'
                pc.dumpJson(filename, scopus_object, filetype='raw/scopus')
                return scopus_object
            except:
                print('No scopus ID found.')
                return None
    except:
        return None
