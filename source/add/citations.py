#! /usr/bin/env python

#Use the elsevier API to get the number of citations a paper has bsed on its PMID.
#Ultimately need to build a GET string like
#http://api.elsevier.com/content/search/scopus?query=PMID(18562177)&apiKey=8024d746590aade6be6856a22a734783&field=citedby-count



def citations(pmids, papers):

    import urllib2
    import json
    import datetime
    import time
    import csv
    import logging

    api_key='8024d746590aade6be6856a22a734783'
    field='citedby-count'
    url='http://api.elsevier.com/content/search/scopus'

    #open the citation cache file
    cached_citations = {}
    try:
        with open('../cache/citations.csv', 'rb') as csvfile:
            f = csv.reader(csvfile)
            for row in f:
                cached_citations[row[0]]={}
                cached_citations[row[0]]['citation_count']=row[1]
                cached_citations[row[0]]['date_downloaded']=row[2]
        csvfile.close()
        logging.info('Citation cache file read in')
    except:
        print 'make file'



    for this_pmid in pmids:

        #read the cache
        try:
            papers[this_pmid]['Extras']['Citations'] = cached_citations[this_pmid]['citation_count']
            logging.info(str(this_pmid)+' in citation cache')
        except:
            #Stick in a small nap so we arent hammering the api too much
            time.sleep(1)

            #try querying with the DOI first - there might not be a DOI
            try:
                request_string=url+'?apiKey='+api_key+'&field=citedby-count&query=DOI('+papers[this_pmid]['doi'][0]+')'
                logging.info(request_string)
                response = urllib2.urlopen(request_string).read()
                t=json.loads(response)

                try:
                    citations = t['search-results']['entry'][0]['citedby-count']
                    #print citations
                    papers[this_pmid]['Extras']['Citations']=citations
                    cached_citations[this_pmid] = {}
                    cached_citations[this_pmid]['citation_count'] = citations
                    cached_citations[this_pmid]['date_downloaded'] = datetime.datetime.now()
                    logging.info('Citation added via DOI')
                except:
                    #there wasnt a number of citations returned, so see if we can catch this.
                    try:
                        error = t['search-results']['entry'][0]['error']
                        if error == 'Result set was empty':
                            logging.info('No citation results from scopus using DOI %s %s',str(papers[this_pmid]['doi']),str(this_pmid))
                    except:
                        #a different error happened!
                        #log this
                        logging.warn('An unexpected error happened getting the citations via DOI!')
                        logging.warn(t)
                        print 'An unexpected error happened getting the citations via DOI!'
                        print request_string
                        print t
                        print t['search-results']['entry'][0]['error']
            except:
                logging.info('No DOI for PMID= '+this_pmid)

            #The above could have failed a couple of points - no DOI or nothing returned from a DOI query
            try:
                papers[this_pmid]['Extras']['Citations']
            except:
                try:
                    #Now try with a PMID
                    request_string=url+'?apiKey='+api_key+'&field=citedby-count&query=PMID('+this_pmid+')'
                    logging.info(request_string)
                    response = urllib2.urlopen(request_string).read()
                    t=json.loads(response)

                    #sometimes this returns multiple entries e.g. 22935244

                    try:
                        citations = t['search-results']['entry'][0]['citedby-count']
                        #print citations
                        papers[this_pmid]['Extras']['Citations']=citations
                        cached_citations[this_pmid] = {}
                        cached_citations[this_pmid]['citation_count'] = citations
                        cached_citations[this_pmid]['date_downloaded'] = datetime.datetime.now()
                        logging.info('Citation added via PMID')
                    except:
                        #there wasnt a number of citations returned, so see if we can catch this.
                        try:
                            error = t['search-results']['entry'][0]['error']
                            if error == 'Result set was empty':
                                #log this
                                logging.info('No citation results from scopus for '+str(this_pmid))
                                #print 'No citations'
                        except:
                            #a different error happened!
                            #log this
                            logging.warn('An unexpected error happened getting the citations!')
                            logging.warn(t)
                            print 'An unexpected error happened getting the citations!'
                            print request_string
                            print t

                except:
                    pass

            try:
                papers[this_pmid]['Extras']['Citations']
            except:
                #If we get here then there is no citation.
                logging.warn('No citations found for %s.',str(this_pmid))

    csvfile = open('../cache/citations.csv', 'wb')
    citation_file =csv.writer(csvfile)
    for this_citation in cached_citations:
        #print str(cached_citations[this_citation])
        temp_citation_count = cached_citations[this_citation]['citation_count']
        temp_date_downloaded = cached_citations[this_citation]['date_downloaded']
        citation_file.writerow([this_citation, str(temp_citation_count), temp_date_downloaded])
