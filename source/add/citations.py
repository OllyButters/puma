#! /usr/bin/env python

#Use the elsevier API to get the number of citations a paper has bsed on its PMID.
#Ultimately need to build a GET string like
#http://api.elsevier.com/content/search/scopus?query=PMID(18562177)&apiKey=8024d746590aade6be6856a22a734783&field=citedby-count



def citations(pmids, papers):
    
    import urllib2
    import json
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
                cached_citations[row[0]]=row[1]
        csvfile.close()
        logging.info('Citation cache file read in')
    except:
        print 'make file'



    for this_pmid in pmids:

        #read the cache
        try:
            papers[this_pmid]['Extras']['Citations'] = cached_citations[this_pmid]
            logging.info(str(this_pmid)+' not in citation cache')
        except:
            #Stick in a small nap so we arent hammering the api too much
            time.sleep(1)
            request_string=url+'?apiKey='+api_key+'&field=citedby-count&query=PMID('+this_pmid+')'
            logging.info(request_string)
            response = urllib2.urlopen(request_string).read()
            t=json.loads(response)

            #sometimes this returns multiple entries e.g. 22935244
            
            try:
                citations = t['search-results']['entry'][0]['citedby-count']
                #print citations
                papers[this_pmid]['Extras']['Citations']=citations
                cached_citations[this_pmid]= citations
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
                    print 'An unexpected error happened getting the citations!')
                    print t

    csvfile = open('../cache/citations.csv', 'wb')
    citation_file =csv.writer(csvfile)
    for this_citation in cached_citations:
        #print str(cached_citations[this_citation])
        temp = cached_citations[this_citation]
        citation_file.writerow([this_citation, str(temp)])


