#! /usr/bin/env python

import csv
import json

############################################################
#Have all the data now, so do something with it

############################################################
#Build a list of all papers by year
def build_yearly(pmids, papers):
    print "\n###HTML yearly list###"

    yearly_papers = {}
    html_file = open('../html/yearly.html', 'w')

    #Build the text needed for each paper
    for this_pmid in pmids:
        try:
            #Paper title as a link
            html='<span style="text-decoration: underline; font-weight:bold;">'+papers[this_pmid]['ArticleTitle']+'</span><br/>'
            
            #Abstract text - probably too long to go on this page
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

            #Make sure there is a dict item for this year
            if this_year not in yearly_papers:
                yearly_papers[this_year] = list()

            temp = yearly_papers[this_year]
            temp.append({this_pmid:html})
            yearly_papers[this_year]=temp
        except:
            pass

    #Output the info into an HTML file
    #For each year dict item
    for this_year in sorted(yearly_papers, key=yearly_papers.get, reverse=True):
        #Check there is some data for this year - not all do
        if len(yearly_papers[this_year])==0:
            continue
        heading='<a name="'+this_year+'"></a>'
        heading+='<h1>'+this_year+'</h1>'
        print >>html_file,heading
        #This is a list
        for this_item in yearly_papers[this_year]:
            temp=this_item.values()
            print >>html_file,temp[0].encode('utf-8')
            

                    

############################################################
#Build a list of all mesh keywords 
def build_mesh(pmids, papers):

    
    print "\n###HTML - mesh###"

    mesh_papers={}
    html_file = open('../html/mesh.html', 'w')

    #Build the text needed for each paper
    for this_pmid in pmids:
        try: 
            for this_mesh in papers[this_pmid]['MeshHeadingList']:          
                if this_mesh['DescriptorName'] not in mesh_papers:
                    mesh_papers[this_mesh['DescriptorName']] = list()
                mesh_papers[this_mesh['DescriptorName']].append(this_pmid)
        except:
            pass

    #print mesh_papers

    #Make a JSON file for each mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers:
        file_name='../html/mesh/'+this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers[this_mesh], indent=4))
        fo.close()
    

    #Make a page with the headings on it
    print >>html_file,'<ul>'
    for this_mesh in sorted(mesh_papers):
        temp = '<li><a href="../html/mesh/'+this_mesh+'">'+this_mesh+'</a></li>'
        print >>html_file,temp
    print >>html_file,'</ul>'


############################################################
#Build a summary page 
def build_mesh(pmids, papers):

    print "\n###HTML - summary###"

    summary={}
    html_file = open('../html/index.html', 'w')

    #Build the text needed for each paper
    for this_pmid in pmids:
        try: 
            this_year=papers[this_pmid]['ArticleDateYear']
            #Make sure there is a dict item for this year
            if this_year not in summary:
                summary[this_year] = {'num_papers':0, 'cumulative':0}

            summary[this_year]['num_papers'] = summary[this_year]['num_papers']+1
        except:
            pass

    #Calculate the cumulative number of papers published
    for this_year in sorted(summary, reverse=False):
        try:
            summary[this_year]['cumulative']=summary[this_year]['num_papers']+summary[str(int(this_year)-1)]['cumulative']
        except:
            summary[this_year]['cumulative']=summary[this_year]['num_papers']

    #print summary

    #Make a page with the headings on it
    print >>html_file,'<table border="1">'
    print >>html_file,'<tr><th>Year</th><th>Number published</th><th>Cumulative</th></tr>'
    for this_year in sorted(summary, reverse=True):
        temp = '<tr><td><a href="yearly.html#'+this_year+'">'+this_year+'</a></td>'
        temp += '<td>'+str(summary[this_year]['num_papers'])+'</td>'
        temp += '<td>'+str(summary[this_year]['cumulative'])+'</td></tr>'
        print >>html_file,temp
    print >>html_file,'</table>'
