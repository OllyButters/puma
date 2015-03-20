#! /usr/bin/env python

#Use the elsevier API to get the number of citations a paper has bsed on its PMID.
#Ultimately need to build a GET string like
#http://api.elsevier.com/content/search/scopus?query=PMID(18562177)&apiKey=8024d746590aade6be6856a22a734783&field=citedby-count

def citations(pmids, papers):
    
    import urllib2
    import json
    
    
    api_key='8024d746590aade6be6856a22a734783'
    PMID='21302344'
    field='citedby-count'
    url='http://api.elsevier.com/content/search/scopus'
    

    for this_pmid in pmids:
        try:
            request_string=url+'?apiKey='+api_key+'&field=citedby-count&query=PMID('+this_pmid+')'
            print request_string
            response = urllib2.urlopen(request_string).read()
            t=json.loads(response)
            
            #print t
            #print '#####################'
            #print t['search-results']['entry'][0]['citedby-count']
            
            citations = t['search-results']['entry'][0]['citedby-count']
            print citations
            papers[this_pmid]['Extras']['Citations']=citations
        except:
            print 'broken'
    

