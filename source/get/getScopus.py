import urllib.request, urllib.error, urllib.parse
import json
import logging

from config import config
from . import papersCache as pc


################################################################################
# Use the elsevier API to get extra info about the paper (including citations).
# Try using the PMID first, if nothing returned then try using the DOI.
################################################################################
def getScopus(zotero_ID, PMID, DOI):

    url = 'https://api.elsevier.com/content/search/scopus'

    try:
        # query scopus with a pmid
        if PMID != "":
            request_string = url + '?apiKey=' + config.scopus_api_key + '&query=PMID(' + str(PMID) + ')'
            logging.info(request_string)
            try:
                response = urllib.request.urlopen(request_string).read()
            except urllib.error.HTTPError as e:  # Catch HTTPError specifically
                if e.code == 429:
                    print("Too many Scopus requests")
                    print(str(e))
                    logging.warning("Too many Scopus requests")
                    logging.warning(str(e))
                    return "QUOTAEXCEEDED"
                
                print("Uncaught HTTP error" + str(e))
                logging.warning("Uncaught HTTP error %s", str(e))
                return
            except Exception as e:  # Catch other exceptions
                print("Error: " + str(e))
                logging.warning("Error: %s", str(e))
                return

            scopus_object = json.loads(response)

            # Check to see if this is just an errorlog
            try:
                error = scopus_object['search-results']['entry'][0]['error']
                if error == 'Result set was empty':
                    logging.info('Result set empty for %s', str(PMID))
                    logging.warning(scopus_object)
                    del(scopus_object)
            except:
                pass

            if scopus_object:
                logging.info('Scopus data got via PMID.')
            else:
                raise Exception('No scopus data returned.')
        else:
            raise Exception('No PMID')
    except:
        print('Unexpected error with PMID.')
        logging.warning('Unexpected error with PMID.')

        # querying with PMID failed, so try DOI.
        try:
            if DOI != "":
                try:
                    request_string = url + '?apiKey=' + config.scopus_api_key + '&query=DOI(' + DOI + ')'
                    logging.info(request_string)

                    try:
                        response = urllib.request.urlopen(request_string).read()
                    except urllib.error.HTTPError as e:  # Catch HTTPError specifically
                        if e.code == 429:
                            print("Too many Scopus requests")
                            print(str(e))
                            logging.warning("Too many Scopus requests")
                            logging.warning(str(e))
                            return "QUOTAEXCEEDED"
                        else:
                            print("Uncaught error" + str(e))
                            logging.warning("Uncaught error" + str(e))
                            return
                    except Exception as e:  # Catch other exceptions
                        print("Error: " + str(e))
                        logging.warning("Error: " + str(e))
                        return

                    scopus_object = json.loads(response)

                    # Check to see if this is just an errorlog
                    try:
                        error = scopus_object['search-results']['entry'][0]['error']
                        if error == 'Result set was empty':
                            logging.info('Result set empty for ' + str(DOI))
                            logging.warning(scopus_object)
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
