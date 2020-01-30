#!/usr/bin/env python3

# import urllib2
import urllib.request, urllib.error, urllib.parse
import json
import datetime
import time
import csv
import logging

import config.config as config


def citations(papers, api_key, citation_max_life, force_update):

    ############################################################################

    # === Europe PMC ===
    print('Doing Europe PMC Citations')
    cached_citations = {}

    # If the force_update flag is True then there is no point reading in the cache
    if force_update is not True:
        logging.info('Reading citation cache in.')
        try:
            # removed 'rb' as no longer need to read bytes (but tbc)
            with open(config.cache_dir + '/citations_europePMC.csv', 'r') as csvfile:
                f = csv.reader(csvfile)
                for row in f:
                    # Parse the date the citation was cached
                    date_downloaded = datetime.datetime.strptime(row[2], "%Y-%m-%d %H:%M:%S.%f")

                    # If the citation is younger than (today - citation_max_life)
                    # then use it. Not using it means we will download it again.
                    if abs(datetime.datetime.now() - date_downloaded) < datetime.timedelta(days=citation_max_life):
                        cached_citations[row[0]] = {}
                        cached_citations[row[0]]['citation_count'] = row[1]
                        cached_citations[row[0]]['date_downloaded'] = row[2]
                        logging.debug('Cache ok for: ' + row[0] + ' ' + row[1] + '' + row[2])
                    else:
                        logging.debug('Cache too old for: ' + row[0] + ' ' + row[1] + '' + row[2])
            csvfile.close()
            logging.info('Citation cache file read in')
        except:
            logging.info('No PMC citation cache file present. Will build a new one.')

    # Get Citation Data
    counter = 1
    for this_paper in papers:
        try:
            this_paper['clean']['citations']['PMC']['count'] = cached_citations[this_paper['IDs']['hash']]['citation_count']
            logging.info(str(this_paper['IDs']['hash'])+" in EuropePMC citation cache (" + str(counter) + "/" + str(len(papers)) + ")")
        except:
            try:
                request_string = 'http://www.ebi.ac.uk/europepmc/webservices/rest/search?format=JSON&query=' + this_paper['IDs']['PMID']
                # response = urllib2.urlopen(request_string).read()
                response = urllib.request.urlopen(request_string).read()
                t = json.loads(response)

                time.sleep(0.4)

                citations = t["resultList"]["result"][0]["citedByCount"]
                this_paper['clean']['citations']['PMC']['count'] = citations
                this_paper['clean']['citations']['PMC']['date_downloaded'] = str(datetime.datetime.now())

                cached_citations[this_paper['IDs']['hash']] = {}
                cached_citations[this_paper['IDs']['hash']]['citation_count'] = citations
                cached_citations[this_paper['IDs']['hash']]['date_downloaded'] = datetime.datetime.now()

                logging.info(str(this_paper['IDs']['hash'])+" fetched EuropePMC citation count (" + str(counter) + "/" + str(len(papers)) + ")")

            except:
                pass

        counter += 1

    # Write to file
    csvfile = open(config.cache_dir + '/citations_europePMC.csv', 'w')
    citation_file = csv.writer(csvfile)
    for this_citation in cached_citations:
        temp_citation_count = cached_citations[this_citation]['citation_count']
        temp_date_downloaded = cached_citations[this_citation]['date_downloaded']
        citation_file.writerow([this_citation, str(temp_citation_count), temp_date_downloaded])
