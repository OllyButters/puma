#! /usr/bin/env python

import json
import sys

# Version 2 of the html pages (New site that matches bristol site)

############################################################
# Have all the data now, so do something with it
############################################################

site_title = "ALSPAC Data Set Publications"

# === Common Page Features ===


def build_common_body(breadcrumb, nav_path, body):
    # nav_path used for changes to relative pathing depending on the page (ie Home does not need anything but next levels need leading ../)
    html = "<body " + body + ">"

    html += "<div class='uob-header-container'>"
    html += "<div class='uob-header width-master' role='banner'>"

    html += "<div class='title_stop'></div>"
    html += "<div id='uoblogo'><a accesskey='1' title='University of Bristol homepage' href='http://www.bristol.ac.uk/''><span>University of Bristol</span></a></div>"
    html += "<div class='maintitle' id='maintitle1'>"
    html += "<span id='title1'><a href='" + nav_path + "index.html'>" + site_title + "</a></span>"
    html += "</div>"
    html += "</div>"
    html += "</div>"

    html += '<div id="uobcms-wrapper" class="width-master">'
    html += '<div id="uobcms-col1" role="navigation">'
    html += '<!--htdig_noindex-->'
    html += '<h4 class="navtitle">'
    html += '<!-- navigation object : navigation title -->'
    html += '<a href="http://www.bristol.ac.uk/alspac/">Avon Longitudinal Study of Parents and Children</a>'
    html += '</h4>'
    html += '<div class="before-navgroup">'
    html += '<!-- navigation object : navigation top -->'

    html += '</div>'
    html += '<!-- navigation object : navigation -->'

    html += '<ul class="navgroup">'
    html += '<li><a href="' + nav_path + 'index.html">Home</a></li>'
    html += '<li><a href="' + nav_path + 'papers/index.html">Papers List</a></li>'
    html += '<li><a href="' + nav_path + 'all_keywords/index.html">All Keywords</a></li>'
    html += '<li><a href="' + nav_path + 'major_keywords/index.html">Major Keywords</a></li>'
    html += '<li><a href="' + nav_path + 'map/index.html">Publication Map</a></li>'
    html += '<li><a href="' + nav_path + 'metrics/index.html">Study Metrics</a></li>'
    html += '</ul>'

    html += '<div class="after-navgroup">'
    html += '<!-- navigation object : navigation bottom -->'
    html += '<!-- start navigation : additional logo -->'
    html += '<div class="logo-additional">'
    html += '<a href="http://www.bristol.ac.uk/alspac/25/"><img src="http://www.bristol.ac.uk/media-library/sites/alspac/images/alspac-25-logo.png" alt="" width="279" height="375" /></a>&zwnj;'
    html += '</div>'
    html += '</div>'
    html += '</div>'

    html += breadcrumb

    html += '<div id="uobcms-content" role="main">'

    return html


def build_common_foot():

    html = '</div>'
    html += '</div>'

    html += '<div class="feedback-container width-master clear"><div class="page-feedback small"></div></div>'
    # html += '<div class="uob-mega-footer-container clearfix"><div class="uob-mega-footer width-master"><div class="uob-mega-footer-col uob-mega-footer-col-1"></div></div></div>'
    html += '<div class="foot clearfix"></div>'
    html += '</body>'
    html += '</html>'

    return html

############################################################
# Home page with summary of years
############################################################


def build_home(papers):

    import shutil

    print "\n###HTML - Home###"

    summary = {}
    html_file = open('../html/index.html', 'w')
    data_file = open('../html/data.js', 'w')

    # Copy CSS files
    shutil.copyfile('html/templates/style_main.css', '../html/css/style_main.css')
    shutil.copyfile('html/templates/uobcms_corporate.css', '../html/css/uobcms_corporate.css')
    shutil.copyfile('html/templates/colour_scheme.css', '../html/css/colour_scheme.css')

    # Put html together for this page
    temp = '<html>'

    # html head
    temp += '<head>'
    temp += '<title>' + site_title + '</title>'
    temp += '<link rel="stylesheet" href="css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="css/style_main.css">'
    temp += '</head>'

    temp += build_common_body("", "", "")
    temp += '<h1 id="pagetitle">Summary by Year</h1>'

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
    # shutil.copyfile('html/templates/plot.html','../html/plot.html')

    # Cohort-Rating calculation
    cr_current_year = 2016.0
    cr_sum = 0.0

    # print summary
    # Make a page with the headings on it
    print >>html_file, '<table>'
    print >>html_file, '<tr><th>Year</th><th>Number published</th><th>Cumulative</th><th>UoB #</th><th>UoB %</th><th>Citations for papers published in this year</th><th>Cumulative citations for papers published in this year</th></tr>'
    for this_year in sorted(summary, reverse=True):
        # Skip the years where nothing was published
        if summary[this_year]['num_papers'] == 0:
            continue

        cr_sum += float(summary[this_year]['citations']) / (cr_current_year - float(this_year))

        # Build the table
        temp = '<tr><td><a href="papers/index.html#'+str(this_year)+'">'+str(this_year)+'</a></td>'
        temp += '<td>'+intWithCommas(summary[this_year]['num_papers'])+'</td>'
        temp += '<td>'+str(summary[this_year]['cumulative'])+'</td>'
        temp += '<td>'+intWithCommas(summary[this_year]['uob'])+'</td>'
        temp += '<td>'+str(int(100*summary[this_year]['uob']/summary[this_year]['num_papers']))+'</td>'
        temp += '<td>'+intWithCommas(summary[this_year]['citations'])+'</td>'
        temp += '<td>'+intWithCommas(summary[this_year]['cumulative_citations'])+'</td></tr>'
        print >>html_file, temp
    print >>html_file, '</table>'

    temp = build_common_foot()
    print >>html_file, temp

    cr_sum = cr_sum / len(papers)
    return cr_sum
############################################################
# Home page with summary of years
############################################################


def build_papers(papers):

    import shutil
    print "\n###HTML papers list###"

    yearly_papers = {}
    html_file = open('../html/papers/index.html', 'w')

    shutil.copyfile('html/templates/altmetric.png', '../html/papers/altmetric.png')

    # Put html together for this page
    temp = '<html>'

    # html head
    temp += '<head>'
    temp += '<title>' + site_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Papers List</p>', "../", "")

    temp += '<h1 id="pagetitle">Papers List</h1>'
    temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

    print >>html_file, temp

    # Build the text needed for each paper
    for this_paper in papers:

        for x in this_paper:
            print (x, ':', this_paper[x])

        try:
            html = ''

            # altmetric data
            try:
                html += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>'
                # html += '<br/><img class="altmetric" src="' + papers[this_pmid]['altmetric_badge'] + '" alt="Badge"><a href="'+papers[this_pmid]['details_url']+'">Altmetric</a>'
            except:
                pass

            # Paper title as a link
            html += '<span style="text-decoration: underline; font-weight:bold;">' + this_paper['title'] + '</span><br/>'

            # Abstract text - probably too long to go on this page
            # html += papers[this_pmid]['AbstractText'][0]+'<br/>'

            # Authors
            authors = []
            for this_author in this_paper['author']:
                # Some author lists have a collective name. Ignore this.
                # Some people don't actually have initials. eg wraight in pmid:18454148
                try:
                    authors.append(this_author['family'] + ', ' + this_author['Initials'])
                except:
                    pass

            html += '; '.join(authors)
            html += '<br/>'

            # Journal volume
            try:
                html += this_paper['Journal'] + ' Vol ' + this_paper['JournalVolume'] + '<br/>'
            except:
                pass

            # PMID
            try:
                html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + str(this_paper['IDs']['PMID'])+'">'+str(this_paper['IDs']['PMID']) + '</a>'
            except:
                pass

            # Zotero
            try:
                html += '&nbsp;Zotero: <a href="' + this_paper['IDs']['zotero'] + '">' + this_paper['IDs']['zotero'] + '</a>'
            except:
                pass

            # DOI
            try:
                html += '&nbsp;DOI: <a href="http://doi.org/' + this_paper['IDs']['DOI'] + '">' + this_paper['IDs']['DOI'] + '</a>'
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

            # Make sure there is a dict item for this year
            if this_year not in yearly_papers:
                yearly_papers[this_year] = list()

            temp = yearly_papers[this_year]
            temp.append({this_paper['IDs']['hash']: html})
            yearly_papers[this_year] = temp
        except:
            print 'Failing on ' + this_paper['IDs']['hash']
            print sys.exc_info()
            pass

    # Output the info into an HTML file
    # For each year dict item
    # for this_year in sorted(yearly_papers, key=yearly_papers.get, reverse=True):
    for this_year in sorted(yearly_papers, reverse=True):
        # Check there is some data for this year - not all do
        if len(yearly_papers[this_year]) == 0:
            continue
        heading = '<a name="' + this_year + '"></a>'
        heading += '<h1>' + this_year + '</h1>'
        print >>html_file, heading
        # This is a list
        for this_item in yearly_papers[this_year]:
            temp = this_item.values()
            print >>html_file, temp[0].encode('utf-8')

    temp = build_common_foot()
    print >>html_file, temp


############################################################
# Build a list of all mesh keywords
############################################################
def build_mesh(papers):
    import os.path

    print "\n###HTML - mesh###"

    mesh_papers_all = {}
    mesh_papers_major = {}
    html_file_all = open('../html/all_keywords/index.html', 'w')
    html_file_major = open('../html/major_keywords/index.html', 'w')

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

    # Print mesh_papers

    # Make a JSON file for each mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers_all:
        file_name = '../html/mesh/all_' + this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers_all[this_mesh], indent=4))
        fo.close()

    # Make a JSON file for each major mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers_major:
        file_name = '../html/mesh/major_'+this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers_major[this_mesh], indent=4))
        fo.close()

    # Put html together for this page
    temp = '<html>'

    # html head
    temp += '<head>'
    temp += '<title>' + site_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; All Keywords</p>', "../", "")

    temp += '<h1 id="pagetitle">All Keywords</h1>'

    print >>html_file_all, temp

    # Make a page with ALL the headings on it
    print >>html_file_all, '<ul>'
    for this_mesh in sorted(mesh_papers_all):
        temp = '<li><a href="../mesh/'+this_mesh+'/index.html">'+this_mesh+'</a></li>'
        print >>html_file_all, temp
    print >>html_file_all, '</ul>'

    temp = build_common_foot()
    print >>html_file_all, temp

    # Put html together for this page
    temp = '<html>'

    # html head
    temp += '<head>'
    temp += '<title>' + site_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Major Keywords</p>', "../", "")

    temp += '<h1 id="pagetitle">Major Keywords</h1>'

    print >>html_file_major, temp

    # Make a page with the MAJOR headings on it
    print >>html_file_major, '<ul>'
    for this_mesh in sorted(mesh_papers_major):
        temp = '<li><a href="../mesh/'+this_mesh+'/index.html">'+this_mesh+'</a></li>'
        print >>html_file_major, temp
    print >>html_file_major, '</ul>'

    temp = build_common_foot()
    print >>html_file_major, temp

    import codecs

    # Make papers list for headings pages
    for this_mesh in mesh_papers_all:

        if (not os.path.exists('../html/mesh/'+this_mesh)):
            os.mkdir('../html/mesh/'+this_mesh)

        file_name = '../html/mesh/' + this_mesh + '/index.html'
        with codecs.open(file_name, 'wb', "utf-8") as fo:

            # Put html together for this page
            temp = '<html>'

            # html head
            temp += '<head>'
            temp += '<title>' + site_title + '</title>'
            temp += '<link rel="stylesheet" href="../../css/uobcms_corporate.css">'
            temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'
            temp += '<link rel="stylesheet" href="../../css/style_main.css">'
            temp += '</head>'

            temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Keyword &gt; ' + this_mesh + '</p>', "../../", "")

            temp += '<h1 id="pagetitle">Keyword - ' + this_mesh + '</h1>'

            print >>fo, temp

            # Build the text needed for each paper
            for this_paper in mesh_papers_all[this_mesh]:

                try:
                    html = ''

                    # altmetric data
                    try:
                        html += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>'
                        # html += '<br/><img class="altmetric" src="' + papers[this_pmid]['altmetric_badge'] + '" alt="Badge"><a href="'+papers[this_pmid]['details_url']+'">Altmetric</a>'
                    except:
                        pass

                    # Paper title as a link
                    html += '<span style="text-decoration: underline; font-weight:bold;">' + this_paper['title'] + '</span><br/>'

                    # Abstract text - probably too long to go on this page
                    # html += papers[this_pmid]['AbstractText'][0]+'<br/>'

                    # Authors
                    authors = []
                    for this_author in this_paper['author']:
                        # Some author lists have a collective name. Ignore this.
                        # Some people don't actually have initials. eg wraight in pmid:18454148
                        try:
                            authors.append(this_author['family'] + ', ' + this_author['Initials'])
                        except:
                            pass

                    html += '; '.join(authors)
                    html += '<br/>'

                    # Journal volume
                    try:
                        html += this_paper['Journal'] + ' Vol ' + this_paper['JournalVolume'] + '<br/>'
                    except:
                        pass

                    # PMID
                    try:
                        html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + str(this_paper['IDs']['PMID'])+'">'+str(this_paper['IDs']['PMID']) + '</a>'
                    except:
                        pass

                    # DOI
                    try:
                        html += '&nbsp;DOI: <a href="http://doi.org/' + this_paper['IDs']['DOI'] + '">' + this_paper['IDs']['DOI'] + '</a>'
                    except:
                        pass

                    # citation count
                    try:
                        html += '&nbsp; Citations: ' + this_paper['Extras']['Citations']
                    except:
                        pass

                    # Add an extra line break at the end
                    html += '<br/><br/>'

                    fo.write(html)
                    print "NOT FAIL"
                except:
                    print 'Failing on ' + this_paper['IDs']['hash']
                    print sys.exc_info()
                    pass

            temp = build_common_foot()
            print >>fo, temp

        fo.close()

###########################################################
# Build a google map based on the lat longs provided before.
###########################################################


def build_google_map(papers):

    import shutil

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
        kml += '["' + this_info['name'] + '",' + this_info['lat'] + ',' + this_info['long'] + '],'
    kml += ']'

    kml_file = open('../html/map/map.kml', 'w')
    print >>kml_file, kml

    html_file = open('../html/map/index.html', 'w')

    # Put html together for this page
    temp = '<html>'

    # html head
    temp += '<head>'
    temp += '<title>' + site_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="http://maps.googleapis.com/maps/api/js?sensor=false"></script>'
    temp += '<script type="text/javascript" src="map.kml"></script>'
    temp += '<script type="text/javascript" src="map.js"></script>'

    shutil.copyfile('html/templates/map.js', '../html/map/map.js')
    shutil.copyfile('html/templates/loading.gif', '../html/map/loading.gif')
    shutil.copyfile('html/templates/map.css', '../html/css/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications Map</p>', "../", "onload='initialize()'")

    temp += '<h1 id="pagetitle">Publications Map</h1>'

    temp += "<div class='loading'><img src='loading.gif'></div>"
    temp += "<div id='map_canvas'></div>"

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp

###########################################################
# Build metrics page
###########################################################


def intWithCommas(x):
    if type(x) not in [type(0), type(0L)]:
        raise TypeError("Parameter must be an integer.")
    if x < 0:
        return '-' + intWithCommas(-x)
    result = ''
    while x >= 1000:
        x, r = divmod(x, 1000)
        result = ",%03d%s" % (r, result)
    return "%d%s" % (x, result)


def build_metrics(papers, cohort_rating):

    import shutil

    html_file = open('../html/metrics/index.html', 'w')

    # Put html together for this page
    temp = '<html>'

    # html head
    temp += '<head>'
    temp += '<title>' + site_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    shutil.copyfile('html/templates/metrics.js', '../html/metrics/metrics.js')

    temp += '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    temp += '<script type="text/javascript" src="../data.js"></script>'
    temp += '<script type="text/javascript" src="../map/map.js"></script>'
    temp += '<script type="text/javascript" src="metrics.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Study Metrics</p>', "../", "")

    temp += '<h1 id="pagetitle">Study Metrics</h1>'

    temp += "<p>{Explanation of study metrics}</p>"

    # Metric calculations
    total_publications = len(papers)
    total_citations = 0
    paper_citations = []
    c20_index = 0

    study_start_year = 1991
    study_current_year = 2016
    study_duration = study_current_year - study_start_year

    c_index_bound = 100

    for this_paper in papers:
        # add the citations for paper to total
        try:

            cit = int(this_paper['Extras']['Citations'])
            total_citations += cit
            paper_citations.append(cit)

            # increment c20-index if more that 20 citations
            if cit >= c_index_bound:
                c20_index += 1
        except:
            pass

    average_citations = float(total_citations)/float(total_publications)
    i20_index_per_year = float(c20_index)/float(study_duration)

    # cal h-index
    paper_citations.sort(reverse=True)
    h_index = 0
    cits_so_far = 0

    for x in range(1, total_publications):
        if x > paper_citations[x]:
            break
        h_index = x
        cits_so_far += paper_citations[x]

    # cal e-index
    # This is probably calculated wrong
    # e_index = math.sqrt(total_citations - cits_so_far)

    # cal g-index
    g_index = 0
    cits_so_far = 0

    for x in range(1, total_publications):
        cits_so_far += paper_citations[x]
        if cits_so_far < x * x:
            break
        g_index = x

    # Ouput Metrics
    temp += "<div class='metric_con'>"
    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>s-index</div>"
    temp += "<div class='metric_value'>" + str(h_index) + "</div>"
    temp += "<div class='metric_description'>s-index is the largest number s such that s publications from a study have at least s citations.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>g-index</div>"
    temp += "<div class='metric_value'>" + str(g_index) + "</div>"
    temp += "<div class='metric_description'>The largest number n of highly cited articles for which the average number of citations is at least n.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Cohort-Rating</div>"
    temp += "<div class='metric_value'>" + str("{0:.3f}".format(round(cohort_rating, 3))) + "</div>"
    temp += "<div class='metric_description'></div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Average Citations Per Publication</div>"
    temp += "<div class='metric_value'>" + str("{0:.2f}".format(round(average_citations, 2))) + "</div>"
    temp += "<div class='metric_description'>The total number of citations divided by the total number of publications.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>c" + str(c_index_bound) + "-index</div>"
    temp += "<div class='metric_value'>" + intWithCommas(c20_index) + "</div>"
    temp += "<div class='metric_description'>The number of publications from a study that have at least " + str(c_index_bound) + " citations.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>c" + str(c_index_bound) + "-index per Study Year</div>"
    temp += "<div class='metric_value'>" + str("{0:.2f}".format(round(i20_index_per_year, 2))) + "</div>"
    temp += "<div class='metric_description'>The c" + str(c_index_bound) + "-index divided by the number of years the study has been running for.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Total Papers</div>"
    temp += "<div class='metric_value'>" + intWithCommas(total_publications) + "</div>"
    temp += "<div class='metric_description'>This is the number of all publications for the study.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Total Citations</div>"
    temp += "<div class='metric_value'>" + intWithCommas(total_citations) + "</div>"
    temp += "<div class='metric_description'>This is the number of citations to all publications for the study.</div>"
    temp += "</div>"
    temp += "</div>"
    temp += "<div class='clear'></div>"

    temp += '<div id="cumulative_div"></div>'
    temp += '<div id="papers_per_year_div"></div>'

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp