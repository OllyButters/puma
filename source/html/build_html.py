#! /usr/bin/env python

import csv
import json
import re
import sys

############################################################
# Have all the data now, so do something with it


############################################################
# Build a list of all papers grouped together by year
############################################################
def build_yearly(papers):
    print "\n###HTML yearly list###"

    yearly_papers = {}
    html_file = open('../html/yearly.html', 'w')

    # Build the text needed for each paper
    for this_paper in papers:
        try:
            # Paper title as a link
            html = '<span style="text-decoration: underline; font-weight:bold;">' + this_paper['title'] + '</span><br/>'

            # Abstract text - probably too long to go on this page
            # html += papers[this_pmid]['AbstractText'][0]+'<br/>'

            # Authors
            authors = []
            for this_author in this_paper['author']:
                # Some author lists have a collective name. Ignore this.
                # Some people don't actually have initials. eg wraight in pmid:18454148
                try:
                    authors.append(this_author['family']+', '+this_author['Initials'])
                except:
                    pass

            html += '; '.join(authors)
            html += '<br/>'

            # Journal volume
            # try:
            # html += papers[0]['Journal']+' Vol '+papers[0]['MedlineCitation']['Article']['Journal']+'<br/>'
            # except:
            #    pass

            # PMID
            try:
                html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/'+str(this_paper['IDs']['PMID'])+'">'+str(this_paper['IDs']['PMID'])+'</a>'
            except:
                pass

            # DOI
            try:
                html += '&nbsp;DOI: <a href="http://doi.org/'+this_paper['IDs']['DOI']+'">'+this_paper['IDs']['DOI']+'</a>'
            except:
                pass

            # citation count
            try:
                html += '&nbsp; Citations: '+this_paper['Extras']['Citations']
            except:
                pass

            # Add an extra line break at the end
            html += '<br/><br/>'

            # Append this paper to the list indexed by the year
            this_year = this_paper['PubmedData']['History'][0]['Year']
            # this_year = this_paper['published-print']['date-parts'][0][0]

            # Make sure there is a dict item for this year
            if this_year not in yearly_papers:
                yearly_papers[this_year] = list()

            temp = yearly_papers[this_year]
            temp.append({this_paper['IDs']['hash']: html})
            yearly_papers[this_year] = temp
        except:
            print 'Failing on '+this_paper['IDs']['hash']
            print sys.exc_info()
            pass

    # Output the info into an HTML file
    # For each year dict item
    # for this_year in sorted(yearly_papers, key=yearly_papers.get, reverse=True):
    for this_year in sorted(yearly_papers, reverse=True):
        # Check there is some data for this year - not all do
        if len(yearly_papers[this_year]) == 0:
            continue
        heading = '<a name="'+str(this_year)+'"></a>'
        heading += '<h1>'+str(this_year)+'</h1>'
        print >> html_file, heading
        # This is a list
        for this_item in yearly_papers[this_year]:
            temp = this_item.values()
            print >> html_file, temp[0].encode('utf-8')


############################################################
# Build a list of all mesh keywords
############################################################
def build_mesh(papers):

    print "\n###HTML - mesh###"

    mesh_papers_all = {}
    mesh_papers_major = {}
    html_file_all = open('../html/mesh_all.html', 'w')
    html_file_major = open('../html/mesh_major.html', 'w')

    # Build a dict of ALL mesh headings with a list of each pmid in each
    for this_paper in papers:
        try:
            # Look at all the mesh headings for this paper
            for this_mesh in this_paper['MedlineCitation']['MeshHeadingList']:
                # If this mesh term is not already in the dict then add it
                if this_mesh['DescriptorName'] not in mesh_papers_all:
                    mesh_papers_all[this_mesh['DescriptorName']] = list()
                mesh_papers_all[this_mesh['DescriptorName']].append(this_paper['IDs']['hash'])
        except:
            pass

    # Build a dict of ONLY MAJOR mesh headings with a list of each pmid in each
    for this_paper in papers:
        try:
            # Look at all the mesh headings for this paper
            for this_mesh in this_paper['MedlineCitation']['MeshHeadingList']:
                # Only interested in majoy topics
                if this_mesh['MajorTopicYN'] == 'Y':
                    # If this mesh term is not in the dict then add it
                    if this_mesh['DescriptorName'] not in mesh_papers_major:
                        mesh_papers_major[this_mesh['DescriptorName']] = list()
                    mesh_papers_major[this_mesh['DescriptorName']].append(this_paper['IDs']['hash'])
        except:
            pass

    # print mesh_papers

    # Make a JSON file for each mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers_all:
        file_name = '../html/mesh/all_'+this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers_all[this_mesh], indent=4))
        fo.close()

    # Make a JSON file for each major mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers_major:
        file_name = '../html/mesh/major_'+this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers_major[this_mesh], indent=4))
        fo.close()

    # Make a page with ALL the headings on it
    print >>html_file_all, '<ul>'
    for this_mesh in sorted(mesh_papers_all):
        temp = '<li><a href="../html/mesh/'+this_mesh+'">'+this_mesh+'</a></li>'
        print >> html_file_all, temp
    print >> html_file_all, '</ul>'

    # Make a page with the MAJOR headings on it
    print >> html_file_major, '<ul>'
    for this_mesh in sorted(mesh_papers_major):
        temp = '<li><a href="../html/mesh/'+this_mesh+'">'+this_mesh+'</a></li>'
        print >> html_file_major, temp
    print >> html_file_major, '</ul>'


############################################################
# Build a summary page
############################################################
def build_summary(papers):

    import shutil

    print "\n###HTML - summary###"

    summary = {}
    html_file = open('../html/index.html', 'w')
    data_file = open('../html/data.js', 'w')

    # Put some links together for this page
    temp = '<a href="yearly.html">All papers</a>&nbsp;'
    temp += '<a href="mesh_all.html">ALL mesh keywords</a>&nbsp;'
    temp += '<a href="mesh_major.html">MAJOR mesh keywords</a>&nbsp;'
    temp += '<a href="map.html">Map</a>&nbsp;<br/>'
    temp += '<a href="plot.html">Plot</a>&nbsp;<br/>'
    print >>html_file, temp

    # Calculate the number of papers for each year
    for this_paper in papers:
        try:
            this_year = this_paper['PubmedData']['History'][0]['Year']
            # Make sure there is a dict item for this year
            if this_year not in summary:
                summary[this_year] = {'num_papers': 0, 'cumulative': 0, 'uob': 0, 'citations': 0, 'cumulative_citations': 0}

            # increment the number of citaitons by one
            summary[this_year]['num_papers'] = summary[this_year]['num_papers']+1

            # add the citations for this paper to the year running total
            try:
                summary[this_year]['citations'] = summary[this_year]['citations'] + int(this_paper['Extras']['Citations'])
            except:
                pass

            # Get number of UoB papers published this year
            if this_paper['Extras']['CleanInstitute'] == 'University of Bristol':
                summary[this_year]['uob'] = summary[this_year]['uob']+1
        except:
            pass

    # Add in some zeros when there is no papers for this year
    years = summary.keys()
    first_year = min(years)
    last_year = max(years)
    for this_year in range(int(first_year), int(last_year)):
        try:
            summary[str(this_year)]['num_papers']
        except:
            summary[str(this_year)] = {'num_papers': 0, 'cumulative': 0, 'uob': 0, 'citations': 0, 'cumulative_citations': 0}

    # Calculate the cumulative number of papers published
    for this_year in sorted(summary, reverse=False):
        try:
            summary[this_year]['cumulative'] = summary[this_year]['num_papers']+summary[str(int(this_year)-1)]['cumulative']
            summary[this_year]['cumulative_citations'] = summary[this_year]['citations']+summary[str(int(this_year)-1)]['cumulative_citations']
        except:
            summary[this_year]['cumulative'] = summary[this_year]['num_papers']
            summary[this_year]['cumulative_citations'] = summary[this_year]['citations']

    ###################################
    # Make a data file that we can plot

    # Cumulative first
    print >>data_file, 'var cumulative =([[\'Year\', \'Number of papers\'],'
    for this_year in sorted(summary, reverse=False):
        print >>data_file, '[\''+this_year+'\','+str(summary[this_year]['cumulative'])+'],'
    print >>data_file, ']);'

    # Number per year now
    print >>data_file, 'var papers_per_year=([[\'Year\', \'Number of papers\'],'
    for this_year in sorted(summary, reverse=False):
        print >>data_file, '[\''+this_year+'\','+str(summary[this_year]['num_papers'])+'],'
    print >>data_file, ']);'

    # Copy the main html page across
    shutil.copyfile('html/templates/plot.html', '../html/plot.html')

    # print summary

    # Make a page with the headings on it
    print >>html_file, '<table border="1">'
    print >>html_file, '<tr><th>Year</th><th>Number published</th><th>Cumulative</th><th>UoB #</th><th>UoB %</th><th>Citations for papers published in this year</th><th>Cumulative citations for papers published in this year</th></tr>'
    for this_year in sorted(summary, reverse=True):
        # Skip the years where nothing was published
        if summary[this_year]['num_papers'] == 0:
            continue

        # Build the table
        temp = '<tr><td><a href="yearly.html#'+str(this_year)+'">'+str(this_year)+'</a></td>'
        temp += '<td>'+str(summary[this_year]['num_papers'])+'</td>'
        temp += '<td>'+str(summary[this_year]['cumulative'])+'</td>'
        temp += '<td>'+str(summary[this_year]['uob'])+'</td>'
        temp += '<td>'+str(int(100*summary[this_year]['uob']/summary[this_year]['num_papers']))+'</td>'
        temp += '<td>'+str(summary[this_year]['citations'])+'</td>'
        temp += '<td>'+str(summary[this_year]['cumulative_citations'])+'</td></tr>'
        print >>html_file, temp
    print >>html_file, '</table>'


###########################################################
# Build a google map based on the lat longs provided before.
###########################################################
def build_google_map(papers):

    import shutil

    # Copy the main html page across
    shutil.copyfile('html/templates/map.html', '../html/map.html')

    info = []
    for this_paper in papers:
        try:
            this_place = {'lat': this_paper['Extras']['LatLong']['lat'], 'long': this_paper['Extras']['LatLong']['long'], 'name': this_paper['Extras']['CleanInstitute']}
            info.append(this_place)
        except:
            pass

    # print info

    kml = "var locations =["
    for this_info in info:
        kml += '["'+this_info['name']+'",'+this_info['lat']+','+this_info['long']+'],'
    kml += ']'

    kml_file = open('../html/map.kml', 'w')
    print >>kml_file, kml


#########################################
# Build a google heat map based
#########################################
def build_google_heat_map():

    print 'Doing heat map'

    data_file = open('../html/map.js', 'w')

    # Open the institute lookup file we have
    geocode = {}
    lat_long = []
    with open('../config/lat_longs.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            try:
                # Check it is not a comment string first.
                if(re.match('#', row[0])):
                    continue

                # If there is a second element in this row then carry on
                lat_long = {'lat': row[1], 'long': row[2]}
                geocode[row[0]] = lat_long
            except:
                pass

    print >>data_file, 'var map =([[\'Latitude\',\'Longitude\',\'Number of papers\',\'Institute\'],'

    # Open the first author institute csv file we made earlier
    with open('../data/first_authors_inst.csv', 'rb') as csvfile:
        f = csv.reader(csvfile)
        for row in f:
            try:
                inst = row[0]
                count = row[1]

                output = '['+geocode[inst]['lat']+','+geocode[inst]['long']+',\''+inst+'\','+count+'],'
                print >>data_file, output
            except:
                pass

    # Make the data file
    print >>data_file, ']);'
