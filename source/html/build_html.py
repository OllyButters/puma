#! /usr/bin/env python

import csv
import json
import sys

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
                #Some author lists have a collective name. Ignore this.
                try: 
                    this_author['CollectiveName']
                    next
                except:
                    #Some people don't actually have initials. eg wraight in pmid:18454148
                    try:
                        authors.append(this_author['LastName']+', '+this_author['Initials'])
                    except:
                        pass

            html += '; '.join(authors)
            html += '<br/>'

            #Journal volume
            try:
                html += papers[this_pmid]['Journal']+' Vol '+papers[this_pmid]['JournalVolume']+'<br/>'
            except:
                pass

            #PMID
            html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/'+str(this_pmid)+'">'+str(this_pmid)+'</a>'

            #DOI
            try:
                html += '&nbsp;DOI: <a href="http://doi.org/'+papers[this_pmid]['doi'][0]+'">'+papers[this_pmid]['doi'][0]+'</a><br/>'
            except:
                pass

            #Add an extra line break at the end
            html += '<br/>'

            #Append this paper to the list indexed by the year
            this_year=papers[this_pmid]['Year']

            #Make sure there is a dict item for this year
            if this_year not in yearly_papers:
                yearly_papers[this_year] = list()


            temp = yearly_papers[this_year]
            temp.append({this_pmid:html})
            yearly_papers[this_year]=temp
        except:
            print 'Failing on '+this_pmid
            print sys.exc_info()
            pass

    #Output the info into an HTML file
    #For each year dict item
    #for this_year in sorted(yearly_papers, key=yearly_papers.get, reverse=True):
    for this_year in sorted(yearly_papers, reverse=True):
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
def build_summary(pmids, papers):

    print "\n###HTML - summary###"

    summary={}
    html_file = open('../html/index.html', 'w')


    #Put some links together for this page
    temp = '<a href="yearly.html">All papers</a>&nbsp;'
    temp += '<a href="mesh.html">Mesh keywords</a>&nbsp;'
    temp += '<a href="map.html">Map</a>&nbsp;<br/>'
    print >>html_file,temp


    #Calculate the number of papers for each year
    for this_pmid in pmids:
        try: 
            this_year=papers[this_pmid]['Year']
            #Make sure there is a dict item for this year
            if this_year not in summary:
                summary[this_year] = {'num_papers':0, 'cumulative':0, 'uob':0}

            summary[this_year]['num_papers'] = summary[this_year]['num_papers']+1

            #Get number of UoB papers published this year
            if papers[this_pmid]['Extras']['CleanInstitute'] == 'University of Bristol':
                summary[this_year]['uob'] = summary[this_year]['uob']+1
        except:
            pass

    #Add in some zeros when there is no papers for this year
    years = summary.keys()
    first_year = min(years)
    last_year = max(years)
    for this_year in range(int(first_year), int(last_year)):
        try:
            summary[str(this_year)]['num_papers']
        except:
            summary[str(this_year)] = {'num_papers':0, 'cumulative':0, 'uob':0}


    #Calculate the cumulative number of papers published
    for this_year in sorted(summary, reverse=False):
        try:
            summary[this_year]['cumulative']=summary[this_year]['num_papers']+summary[str(int(this_year)-1)]['cumulative']
        except:
            summary[this_year]['cumulative']=summary[this_year]['num_papers']




    #print summary

    #Make a page with the headings on it
    print >>html_file,'<table border="1">'
    print >>html_file,'<tr><th>Year</th><th>Number published</th><th>Cumulative</th><th>UoB #</th><th>UoB %</th></tr>'
    for this_year in sorted(summary, reverse=True):
        #Skip the years where nothing was published
        if summary[this_year]['num_papers'] == 0:
            continue
        
        #Build the table
        temp = '<tr><td><a href="yearly.html#'+str(this_year)+'">'+str(this_year)+'</a></td>'
        temp += '<td>'+str(summary[this_year]['num_papers'])+'</td>'
        temp += '<td>'+str(summary[this_year]['cumulative'])+'</td>'
        temp += '<td>'+str(summary[this_year]['uob'])+'</td>'
        temp += '<td>'+str(int(100*summary[this_year]['uob']/summary[this_year]['num_papers']))+'</td></tr>'
        print >>html_file,temp
    print >>html_file,'</table>'



#Build a google map based on the lat longs provided before.
def build_google_map(pmids, papers):

    import shutil

    #Copy the main html page across
    shutil.copyfile('html/templates/map.html','../html/map.html')

    info = []
    for this_pmid in pmids:
        try:
            this_place = {'lat': papers[this_pmid]['Extras']['LatLong']['lat'], 'long': papers[this_pmid]['Extras']['LatLong']['long'], 'name':papers[this_pmid]['Extras']['CleanInstitute']}
            info.append(this_place)
        except:
            pass
        
        
    #print info

    kml = "var locations =["
    for this_info in info:
        kml+= '["'+this_info['name']+'",'+this_info['lat']+','+this_info['long']+'],'
    kml += ']'
    
    kml_file = open('../html/map.kml', 'w')
    print >>kml_file,kml