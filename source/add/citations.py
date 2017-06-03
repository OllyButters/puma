#! /usr/bin/env python

import urllib2
import json
import datetime
import time
import csv
import logging

import config.config as config


# Use the elsevier API to get the number of citations a paper has.
# Try using the PMID first, if nothing returned then try using the DOI.
# Ultimately need to build a GET string like
# http://api.elsevier.com/content/search/scopus?query=PMID(18562177)&apiKey=8024d746590aade6be6856a22a734783&field=citedby-count
# After scopus try pubmedcentral for a citation count.
# There does seem to be a fair bit of repatition between the DOI and PMID code,
# this is because the DOI returned data sometimes has multiple results, so needs
# a bit of extra parsing.
# Note: eid is the ID scopus assigns to each paper it knows about, this is useful
#     for linking to it on the HTML pages later on.
def citations(papers, api_key, citation_max_life, force_update, error_log):

    url = 'http://api.elsevier.com/content/search/scopus'

    # === Scopus ===
    print 'Doing Scopus Citations'

    # open the citation cache file
    cached_citations = {}

    # If the force_update flag is True then there is no point reading in the cache
    logging.info('scopus_force_citation_update is ' + str(force_update))
    if force_update is not True:
        logging.info('Reading citation cache in.')
        try:
            with open(config.cache_dir + '/citations_scopus.csv', 'rb') as csvfile:
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
                        cached_citations[row[0]]['eid'] = row[3]
                        logging.debug('Cache ok for: ' + row[0] + ' ' + row[1] + ' ' + row[2])
                    else:
                        logging.debug('Cache too old for: ' + row[0] + ' ' + row[1] + '' + row[2])
            csvfile.close()
            logging.info('Citation cache file read in')
        except:
            logging.info('No scopus citation cache file present. Will build a new one.')

    number_papers_to_process = len(papers)
    counter = 0
    for this_paper in papers:
        counter = counter + 1
        logging.info('\non # ' + str(counter) + ' of ' + str(number_papers_to_process))

        # read the cache
        try:
            this_paper['clean']['citations']['scopus']['count'] = cached_citations[this_paper['IDs']['hash']]['citation_count']
            this_paper['clean']['citations']['scopus']['date_downloaded'] = cached_citations[this_paper['IDs']['hash']]['date_downloaded']
            try:
                this_paper['IDs']['scopus'] = cached_citations[this_paper['IDs']['hash']]['eid']
            except:
                pass
            logging.info(str(this_paper['IDs']['hash']) + ' in citation cache')
        except:
            # Stick in a small nap so we arent hammering the api too much
            time.sleep(1)

            # Handle Max Quota Reached
            error_number = 0

            # query scopus with a pmid
            if this_paper['IDs']['PMID'] != "":
                try:
                    # Now try with a PMID
                    # request_string = url + '?apiKey=' + api_key + '&field=citedby-count&query=PMID(' + this_paper['IDs']['PMID'] + ')'
                    request_string = url + '?apiKey=' + api_key + '&query=PMID(' + this_paper['IDs']['PMID'] + ')'
                    logging.info(request_string)
                    response = urllib2.urlopen(request_string).read()
                    t = json.loads(response)

                    # sometimes this returns multiple entries e.g. 22935244
                    try:
                        if len(t['search-results']['entry']) > 1:
                            error_log.logErrorPaper("Multiple different citaton counts found for PMID", this_paper)
                        citations = t['search-results']['entry'][0]['citedby-count']
                        this_paper['clean']['citations']['scopus']['count'] = citations
                        this_paper['clean']['citations']['scopus']['date_downloaded'] = datetime.datetime.now()

                        if len(t['search-results']['entry']) == 1:  # Do not cache if multiple results returned
                            cached_citations[this_paper['IDs']['hash']] = {}
                            cached_citations[this_paper['IDs']['hash']]['citation_count'] = citations
                            cached_citations[this_paper['IDs']['hash']]['date_downloaded'] = datetime.datetime.now()
                            try:
                                cached_citations[this_paper['IDs']['hash']]['eid'] = t['search-results']['entry'][0]['eid']
                                this_paper['IDs']['scopus'] = t['search-results']['entry'][0]['eid']
                            except:
                                pass
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
                this_paper['clean']['citations']['scopus']['count']
            except:
                # If we get here then there is no citation.
                logging.warn('No citations found for %s.', str(this_paper['IDs']['hash']))
            # ==================================================

            # ==================================================
            # The above could have failed a couple of points - no PMID or nothing returned from a PMID query
            try:
                # try querying with the DOI first - there might not be a DOI
                if 'count' not in this_paper['clean']['citations']['scopus'] and this_paper['IDs']['DOI'] != "":
                    # request_string = url+'?apiKey='+api_key+'&field=citedby-count&query=DOI('+this_paper['IDs']['DOI']+')'
                    request_string = url + '?apiKey=' + api_key + '&query=DOI(' + this_paper['IDs']['DOI'] + ')'
                    logging.info(request_string)
                    try:
                        response = urllib2.urlopen(request_string).read()
                        t = json.loads(response)
                    except:
                        logging.error('The citation query failed - maybe it timed out?')
                        print 'The citation query failed - maybe it timed out?'
                        # exit()

                    try:
                        if len(t['search-results']['entry']) > 1:
                            error_log.logWarningPaper("Multiple different citation counts found for DOI", this_paper)
                            # Add up citations for DOIs that exactly match
                            citations = 0
                            title = ""
                            # Loops through all returned results looking for matching DOIs
                            for n in range(0, len(t['search-results']['entry'])):
                                if t['search-results']['entry'][n]['prism:doi'] == this_paper['IDs']['DOI'] and (title == "" or title == t['search-results']['entry'][n]['dc:title']):
                                    citations = citations + int(t['search-results']['entry'][n]['citedby-count'])
                                    this_paper['IDs']['scopus'] = t['search-results']['entry'][n]['eid']
                                    title = t['search-results']['entry'][n]['dc:title']
                            this_paper['clean']['citations']['scopus']['count'] = citations

                        elif len(t['search-results']['entry']) == 1:
                            citations = t['search-results']['entry'][0]['citedby-count']
                            this_paper['clean']['citations']['scopus']['count'] = citations
                            this_paper['clean']['citations']['scopus']['date_downloaded'] = datetime.datetime.now()
                            scopus_id = t['search-results']['entry'][0]['eid']
                            this_paper['IDs']['scopus'] = scopus_id

                            # Do not cache if multiple results returned
                            cached_citations[this_paper['IDs']['hash']] = {}
                            cached_citations[this_paper['IDs']['hash']]['citation_count'] = citations
                            cached_citations[this_paper['IDs']['hash']]['date_downloaded'] = datetime.datetime.now()
                            cached_citations[this_paper['IDs']['hash']]['eid'] = scopus_id

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
                    if this_paper['IDs']['DOI'] == "":
                        logging.info('No DOI for = '+this_paper['IDs']['hash'])

                error_number = 0
            except:
                error_number += 1
                print 'An unexpected error happened getting the citations via DOI!'
                if error_number > 10:
                    print 'Check if reached Scopus MAX_QUOTA'

            # ==================================================

    csvfile = open(config.cache_dir + '/citations_scopus.csv', 'wb')
    citation_file = csv.writer(csvfile)
    for this_citation in cached_citations:
        temp_citation_count = cached_citations[this_citation]['citation_count']
        temp_date_downloaded = cached_citations[this_citation]['date_downloaded']
        try:
            temp_eid = cached_citations[this_citation]['eid']
        except:
            temp_eid = ""
        citation_file.writerow([this_citation, str(temp_citation_count), temp_date_downloaded, temp_eid])

    ############################################################################

    # === Europe PMC ===
    print 'Doing Europe PMC Citations'
    cached_citations = {}

    # If the force_update flag is True then there is no point reading in the cache
    if force_update is not True:
        logging.info('Reading citation cache in.')
        try:
            with open(config.cache_dir + '/citations_europePMC.csv', 'rb') as csvfile:
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
                response = urllib2.urlopen(request_string).read()
                t = json.loads(response)

                time.sleep(0.4)

                citations = t["resultList"]["result"][0]["citedByCount"]
                this_paper['clean']['citations']['PMC']['count'] = citations
                this_paper['clean']['citations']['PMC']['date_downloaded'] = datetime.datetime.now()

                cached_citations[this_paper['IDs']['hash']] = {}
                cached_citations[this_paper['IDs']['hash']]['citation_count'] = citations
                cached_citations[this_paper['IDs']['hash']]['date_downloaded'] = datetime.datetime.now()

                logging.info(str(this_paper['IDs']['hash'])+" fetched EuropePMC citation count (" + str(counter) + "/" + str(len(papers)) + ")")

            except:
                # print "No Europe PMC citations count for " + this_paper['IDs']['hash'] + " (" + str(counter) + "/" + str(len(papers)) + ")"
                pass

        counter += 1

    # Write to file
    csvfile = open(config.cache_dir + '/citations_europePMC.csv', 'wb')
    citation_file = csv.writer(csvfile)
    for this_citation in cached_citations:
        temp_citation_count = cached_citations[this_citation]['citation_count']
        temp_date_downloaded = cached_citations[this_citation]['date_downloaded']
        citation_file.writerow([this_citation, str(temp_citation_count), temp_date_downloaded])
