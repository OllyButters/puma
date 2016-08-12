#! /usr/bin/env python

import urllib2
import json
import datetime
import time
import csv
import logging
import sys

import config.config as config


# Use the elsevier API to get the number of citations a paper has bsed on its PMID.
# Ultimately need to build a GET string like
# http://api.elsevier.com/content/search/scopus?query=PMID(18562177)&apiKey=8024d746590aade6be6856a22a734783&field=citedby-count
def citations(papers, api_key, citation_max_life, force_update):

    url = 'http://api.elsevier.com/content/search/scopus'

    print 'Doing citations'

    # open the citation cache file
    cached_citations = {}

    # If the force_update flag is True then there is no point reading in the cache
    logging.info('scopus_force_citation_update is ' + str(force_update))
    if force_update is not True:
        logging.info('Reading citation cache in.')
        try:
            with open(config.cache_dir + '/citations.csv', 'rb') as csvfile:
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
            print("Unexpected error:", sys.exc_info()[0])
            print 'make file'

    number_papers_to_process = len(papers)
    counter = 0
    for this_paper in papers:
        counter = counter + 1
        logging.info('on # '+str(counter)+' of '+str(number_papers_to_process))

        print this_paper['IDs']
        # exit()

        # read the cache
        try:
            this_paper['Extras']['Citations'] = cached_citations[this_paper['IDs']['hash']]['citation_count']
            logging.info(str(this_paper['IDs']['hash'])+' in citation cache')
        except:
            # Stick in a small nap so we arent hammering the api too much
            time.sleep(1)

            # Handle Max Quota Reached
            error_number = 0
            try:
                # try querying with the DOI first - there might not be a DOI
                if this_paper['IDs']['DOI'] != "":
                    request_string = url+'?apiKey='+api_key+'&field=citedby-count&query=DOI('+this_paper['IDs']['DOI']+')'
                    logging.info(request_string)
                    try:
                        response = urllib2.urlopen(request_string).read()
                        print "DOI: " + str(response)
                        t = json.loads(response)
                    except:
                        logging.error('The citation query failed - maybe it timed out?')
                        print 'The citation query failed - maybe it timed out?'
                        # exit()

                    try:
                        citations = t['search-results']['entry'][0]['citedby-count']
                        this_paper['Extras']['Citations'] = citations
                        cached_citations[this_paper['IDs']['hash']] = {}
                        cached_citations[this_paper['IDs']['hash']]['citation_count'] = citations
                        cached_citations[this_paper['IDs']['hash']]['date_downloaded'] = datetime.datetime.now()
                        logging.info('Citation added via DOI')
                    except:
                        # there wasnt a number of citations returned, so see if we can catch this.
                        try:
                            error = t['search-results']['entry'][0]['error']
                            if error == 'Result set was empty':
                                logging.info('No citation results from scopus using DOI %s %s', str(this_paper['IDs']['DOI']), str(this_paper))
                        except:
                            # a different error happened!
                            # log this
                            logging.warn('An unexpected error happened getting the citations via DOI!')
                            logging.warn(t)
                            print 'An unexpected error happened getting the citations via DOI!'
                            print request_string
                            print t
                            print t['search-results']['entry'][0]['error']
                else:
                    logging.info('No DOI for = '+this_paper['IDs']['hash'])

                error_number = 0
            except:
                error_number += 1
                print 'An unexpected error happened getting the citations via DOI!'
                if error_number > 10:
                    print 'Check if reached Scopus MAX_QUOTA'

            # shoud wrap the above up as a fn and run it with doi and pmid separately
            # The above could have failed a couple of points - no DOI or nothing returned from a DOI query
            if 'Citations' not in this_paper['Extras'] and this_paper['IDs']['PMID'] != "":
                try:
                    # Now try with a PMID
                    request_string = url + '?apiKey=' + api_key + '&field=citedby-count&query=PMID(' + this_paper['IDs']['PMID'] + ')'
                    logging.info(request_string)
                    response = urllib2.urlopen(request_string).read()
                    print "PMID: " + str(response)
                    t = json.loads(response)

                    # sometimes this returns multiple entries e.g. 22935244
                    try:
                        citations = t['search-results']['entry'][0]['citedby-count']
                        # print citations
                        this_paper['Extras']['Citations'] = citations
                        cached_citations[this_paper['IDs']['hash']] = {}
                        cached_citations[this_paper['IDs']['hash']]['citation_count'] = citations
                        cached_citations[this_paper['IDs']['hash']]['date_downloaded'] = datetime.datetime.now()
                        logging.info('Citation added via PMID')
                    except:
                        # there wasnt a number of citations returned, so see if we can catch this.
                        try:
                            error = t['search-results']['entry'][0]['error']
                            if error == 'Result set was empty':
                                # log this
                                logging.info('No citation results from scopus for ' + str(this_paper['IDs']['PMID']))
                                # print 'No citations'
                        except:
                            # a different error happened!
                            # log this
                            logging.warn('An unexpected error happened getting the citations!')
                            logging.warn(t)
                            print 'An unexpected error happened getting the citations!'
                            print request_string
                            print t
                except:
                    print 'An unexpected error happened getting the citations via PMID!'

            try:
                this_paper['Extras']['Citations']
            except:
                # If we get here then there is no citation.
                logging.warn('No citations found for %s.', str(this_paper['IDs']['hash']))

    csvfile = open(config.cache_dir + '/citations.csv', 'wb')
    citation_file = csv.writer(csvfile)
    for this_citation in cached_citations:
        temp_citation_count = cached_citations[this_citation]['citation_count']
        temp_date_downloaded = cached_citations[this_citation]['date_downloaded']
        citation_file.writerow([this_citation, str(temp_citation_count), temp_date_downloaded])
