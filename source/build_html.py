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

    sorted_authors = {'2015':[],'2014':[], '2013':[], '2012':[], '2011':[], '2010':[], '2009':[], '2008':[], '2007':[], '2006':[], '2005':[], '2004':[]}
    html_file = open('../html/list.html', 'w')
    for this_pmid in pmids:
        try:
            #Build the text needed for each paper
            #Paper title as a link
            #html='<a href="http://www.ncbi.nlm.nih.gov/pubmed/'+str(this_pmid)+'">'+papers[this_pmid]['ArticleTitle']+'</a><br/>'
            html='<span style="text-decoration: underline; font-weight:bold;">'+papers[this_pmid]['ArticleTitle']+'</span><br/>'
            
            #Abstract text
            #html += papers[this_pmid]['AbstractText'][0]+'<br/>'
            
            #Authors
            authors = []
            for this_author in papers[this_pmid]['AuthorList']:
                authors.append(this_author['LastName']+', '+this_author['Initials'])
                
            html += '; '.join(authors)
            html += '<br/>'

            #Journal volume
            html += papers[this_pmid]['Journal']+' Vol '+papers[this_pmid]['JournalVolume']+'<br/>'

            #PMID
            html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/'+str(this_pmid)+'">'+str(this_pmid)+'</a>'

            #DOI
            html += '&nbsp;DOI: <a href="http://doi.org/'+papers[this_pmid]['doi'][0]+'">'+papers[this_pmid]['doi'][0]+'</a><br/>'

            #Add an extra line break at the end
            html += '<br/>'

            #Append this paper to the list indexed by the year
            this_year=papers[this_pmid]['ArticleDateYear']
            temp = sorted_authors[this_year]
            temp.append({this_pmid:html})
            sorted_authors[this_year]=temp
        except:
            pass
    
#    print sorted_authors

    #For each year dict item
#    for this_year in sorted_authors:
    for this_year in sorted(sorted_authors, key=sorted_authors.get, reverse=True):
        #print >>html_file,len(sorted_authors[this_year])
        if len(sorted_authors[this_year])==0:
            continue
        heading='<h1>'+this_year+'</h1>'
        print >>html_file,heading
        #This is a list
        for this_item in sorted_authors[this_year]:
 #This is a key value pair
            temp=this_item.values()
            #print temp[0]
            print >>html_file,temp[0].encode('utf-8')
            

                    

