#! /usr/bin/env python

import csv

############################################################
#Have all the data now, so do something with it

############################################################
#Build a list of all journals and count frequencies of each
def build_html(pmids, papers):
    print "\n###HTML###"
        
    #print a list of sorted frequencies
#    html_file = open('../html/list.html', 'w')
#    for this_pmid in pmids:
#        html='<a href="http://www.ncbi.nlm.nih.gov/pubmed/'+str(this_pmid)+'">'+papers[this_pmid]['ArticleTitle']+'</a><br/>'
#        for this_author in papers[this_pmid]['AuthorList']:
#        #There are some entries in the author list that are not actually authors e.g. 21379325 last author
 #           try:
 #               html+=this_author['LastName']
 ##           except:
 #               pass
 #       html+='<br/>'
  #      print >>html_file,html.encode('utf-8')

    sorted_authors = {'2014':[], '2013':[]}
    html_file = open('../html/list.html', 'w')
    for this_pmid in pmids:
        try:
            html='<a href="http://www.ncbi.nlm.nih.gov/pubmed/'+str(this_pmid)+'">'+papers[this_pmid]['ArticleTitle']+'</a><br/>'
            this_year=papers[this_pmid]['ArticleDateYear']
            temp = sorted_authors[this_year]
            temp.append({this_pmid:html})
            sorted_authors[this_year]=temp
        except:
            pass
    
    print sorted_authors

    #For each year dict item
    for this_year in sorted_authors:
        heading='<h1>'+this_year+'</h1>'
        print >>html_file,heading
        #This is a list
        for this_item in sorted_authors[this_year]:
            #This is a key value pair
            temp=this_item.values()
            print temp[0]
            print >>html_file,temp[0].encode('utf-8')
            

                    

