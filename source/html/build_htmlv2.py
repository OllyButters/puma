#! /usr/bin/env python

import json
import sys
import shutil
import os.path
import csv
import htmlentities
import time
import datetime
import math
import codecs
import os


import config.config as config

# Version 2 of the html pages (New site that matches bristol's ALSPAC site)

############################################################
# Have all the data now, so do something with it
############################################################

site_second_title = " Data Set Publications"


# === Common Page Features ===
def build_common_body(breadcrumb, nav_path, body):
    # nav_path used for changes to relative pathing depending on the page (ie Home does not need anything but next level down needs leading ../)
    html = "<body " + body + ">"

    html += "<div class='uob-header-container'>"
    html += "<div class='uob-header width-master' role='banner'>"

    html += "<div class='title_stop'></div>"
    html += "<div id='uoblogo'><a accesskey='1' title='" + config.project_details['header_institution'] + "' href='" + config.project_details['header_institution_url'] + "'><span>" + config.project_details['header_institution'] + "</span></a></div>"
    html += "<div class='maintitle' id='maintitle1'>"
    html += "<span id='title1'><a href='" + nav_path + "index.html'>" + config.project_details['name'] + " - " + site_second_title + "</a></span>"
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
    html += '<li><a href="' + nav_path + 'help/index.html">Information</a></li>'
    # html += '<li><a href="' + nav_path + 'papers/index.html">Papers List</a></li>'
    html += '<li><a href="' + nav_path + 'all_keywords/index.html">All Keywords</a></li>'
    html += '<li><a href="' + nav_path + 'major_keywords/index.html">Major Keywords (MeSH)</a></li>'

    html += '<li><a>Maps</a>'
    html += '<ul class="multilevel-linkul-0">'
    html += '    <li><a href="' + nav_path + 'map/index.html">Institutions Map</a></li>'
    html += '    <li><a href="' + nav_path + 'country/index.html">Publications by Country</a></li>'
    html += '    <li><a href="' + nav_path + 'city/index.html">Publications by UK City</a></li>'
    html += '</ul></li>'

    html += '<li><a href="' + nav_path + 'authornetwork/index.html">Author Network</a></li>'
    html += '<li><a href="' + nav_path + 'metrics/index.html">Study Metrics</a></li>'
    html += '<li><a href="' + nav_path + 'wordcloud/index.html">Major Keyword Cloud</a></li>'
    html += '<li><a href="' + nav_path + 'abstractwordcloud/index.html">Abstract Word Cloud</a></li>'
    html += '<li id="error_page_li" style="display:none;"><a href="' + nav_path + 'errorlog/index.html">Error Log</a></li>'
    html += '</ul>'

    # Cookie errorlog display
    html += '<script>function getCookie(cname) {var name = cname + "=";var ca = document.cookie.split(";");for(var i = 0; i <ca.length; i++) { var c = ca[i]; while (c.charAt(0)==" ") {c = c.substring(1);}if (c.indexOf(name) == 0) {return c.substring(name.length,c.length);}}return "";}</script>'
    html += '<script>if ( getCookie("show_error_page") != "" ){ document.getElementById("error_page_li").style.display = "block";}</script>'

    html += '<div class="after-navgroup">'
    html += '<!-- navigation object : navigation bottom -->'
    html += '<!-- start navigation : additional logo -->'
    html += '<div class="logo-additional">'
    html += '<a href="' + config.project_details['side_image_link'] + '"><img src="' + config.project_details['side_image_url'] + '" alt="" width="279" height="375" /></a>&zwnj;'
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
    html += '<div class="foot clearfix"></div>'
    html += '</body>'
    html += '</html>'

    return html


############################################################
# Home page with summary of years
############################################################
def build_home(papers, error_log):

    print "\n###HTML - Home###"

    summary = {}
    missing_year = {'num_papers': 0, 'uob': 0, 'citations': 0}

    html_file = open(config.html_dir + '/index.html', 'w')
    data_file = open(config.html_dir + '/data.js', 'w')

    # Copy CSS files
    shutil.copyfile(config.template_dir + '/style_main.css', config.html_dir + '/css/style_main.css')
    shutil.copyfile(config.template_dir + '/uobcms_corporate.css', config.html_dir + '/css/uobcms_corporate.css')
    shutil.copyfile(config.template_dir + '/colour_scheme.css', config.html_dir + '/css/colour_scheme.css')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="css/style_main.css">'
    temp += '<link rel="stylesheet" href="css/colour_scheme.css">'
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
            summary[this_year]['num_papers'] += 1

            # add the citations for this paper to the year running total
            try:
                summary[this_year]['citations'] += int(this_paper['Extras']['Citations'])
            except:
                pass

            # Get number of UoB papers published this year
            try:
                if this_paper['Extras']['CleanInstitute'] == 'University of Bristol':
                    summary[this_year]['uob'] += 1
            except:
                pass

        except:
            try:
                this_paper['PubmedData']['History'][0]['Year']
            except:
                error_log.logErrorPaper("Date Missing for " + this_paper['IDs']['hash'], this_paper)
                missing_year['num_papers'] += 1
                try:
                    missing_year['citations'] += int(this_paper['Extras']['Citations'])
                except:
                    pass
                    try:
                        if this_paper['Extras']['CleanInstitute'] == 'University of Bristol':
                            missing_year['uob'] += 1
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
        print >>data_file, '[\''+this_year+'\','+str(summary[this_year]['cumulative']) + '],'
    print >>data_file, ']);'

    # Number per year now
    print >>data_file, 'var papers_per_year=([[\'Year\', \'Number of papers\'],'
    for this_year in sorted(summary, reverse=False):
        print >>data_file, '[\''+this_year+'\','+str(summary[this_year]['num_papers']) + '],'
    print >>data_file, ']);'

    # Copy the main html page across
    # shutil.copyfile('html/templates/plot.html','../html/plot.html')

    # Cohort-Rating calculation
    cr_current_year = float(config.metrics_study_current_year)
    print cr_current_year
    cr_sum = 0.0
    cr_data_from = 0

    # print summary
    # Make a page with the headings on it
    print >>html_file, '<table>'
    print >>html_file, '<tr><th>Year</th><th>Number published</th><th>Cumulative</th><th>UoB #</th><th>UoB %</th><th>Citations for papers published in this year</th><th>Cumulative citations for papers published in this year</th></tr>'
    for this_year in sorted(summary, reverse=True):
        # Skip the years where nothing was published
        if summary[this_year]['num_papers'] == 0:
            continue

        cr_year = float(cr_current_year - float(this_year))
        if cr_year < 1:
            cr_year = 1.0

        cr_sum += float(summary[this_year]['citations']) / cr_year
        cr_data_from += summary[this_year]['num_papers']

        # Build the table
        temp = '<tr><td><a href="papers/' + this_year + '/index.html">' + str(this_year) + '</a></td>'
        temp += '<td>' + intWithCommas(summary[this_year]['num_papers']) + '</td>'
        temp += '<td>' + str(summary[this_year]['cumulative']) + '</td>'
        temp += '<td>' + intWithCommas(summary[this_year]['uob']) + '</td>'
        temp += '<td>' + str(int(100*summary[this_year]['uob']/summary[this_year]['num_papers'])) + '</td>'
        temp += '<td>' + intWithCommas(summary[this_year]['citations']) + '</td>'
        temp += '<td>' + intWithCommas(summary[this_year]['cumulative_citations']) + '</td></tr>'
        print >>html_file, temp

        temp = '<tr>'
        temp += '<td style="font-size:12px;font-weight:bold;">UNKNOWN</td>'
        temp += '<td>' + str(missing_year['num_papers']) + '</td>'
        temp += '<td>-</td>'
        temp += '<td>' + str(missing_year['uob']) + '</td>'
        temp += '<td>' + str(int(100*missing_year['uob']/missing_year['num_papers'])) + '</td>'
        temp += '<td>' + str(missing_year['citations']) + '</td>'
        temp += '<td>-</td>'
        temp += '</tr>'
    print >>html_file, temp
    print >>html_file, '</table>'

    temp = "<p>Known publication year for " + intWithCommas(cr_data_from) + " of " + intWithCommas(len(papers)) + " publications</p>"

    temp += build_common_foot()
    print >>html_file, temp

    cr_sum = cr_sum / len(papers)
    return cr_sum, cr_data_from


############################################################
# Home page with summary of years
############################################################
def build_papers(papers):

    print "\n###HTML papers list###"

    yearly_papers = {}
    html_file = open(config.html_dir + '/papers/index.html', 'w')

    shutil.copyfile(config.template_dir + '/altmetric.png', config.html_dir + '/papers/altmetric.png')
    shutil.copyfile(config.template_dir + '/yellow-flag-th.png', config.html_dir + '/papers/yellow-flag-th.png')

    # Read in exec csv
    exec_list = []
    f = open(config.config_dir + "/" + config.project_details['short_name'] + '_exec_members.csv', 'rt')
    try:
        reader = csv.reader(f)
        n = 0
        for row in reader:
            if n > 0:
                exec_list.append(row)
            n += 1
    finally:
        f.close()

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Papers List</p>', "../", "")

    temp += '<h1 id="pagetitle">Papers List</h1>'
    temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

    print >>html_file, temp
    main = "<p>"

    # Build the text needed for each paper
    for this_paper in papers:

        try:
            html = ''

            # altmetric data
            try:
                if this_paper['IDs']['DOI']:
                    html += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>'
            except:
                pass

            # Paper title as a link
            html += '<span style="text-decoration: underline; font-weight:bold;">' + this_paper['title'] + '</span><br/>'

            # Authors
            authors = []
            author_on_exec = False
            for this_author in this_paper['author']:
                # Some author lists have a collective name. Ignore this.
                # Some people don't actually have initials. eg wraight in pmid:18454148
                try:
                    # Check if an author was on the exec Comittee
                    for x in exec_list:
                        # Check if authors name matches
                        if x[2] == this_author['clean']:
                            # Get start and end date of exec membership
                            exec_start = time.mktime(datetime.datetime.strptime(x[0], "%d/%m/%Y").timetuple())
                            exec_end = int(time.time())
                            if not x[1] == "":
                                exec_end = time.mktime(datetime.datetime.strptime(x[1], "%d/%m/%Y").timetuple())

                            # Convert issued date into a timestamp
                            issued_date = this_paper['issued']['date-parts']
                            date = str(issued_date[0][2]) + "/" + str(issued_date[0][1]) + "/" + str(issued_date[0][0])
                            issued_timestamp = time.mktime(datetime.datetime.strptime(date, "%d/%m/%Y").timetuple())

                            # If publication is issued between exec_start and exec_end then flag
                            if issued_timestamp > exec_start and issued_timestamp < exec_end:
                                author_on_exec = True
                                break

                    authors.append(this_author['family'] + ', ' + this_author['given'])
                except:
                    pass

            html += htmlentities.encode('; '.join(authors))
            html += '<br/>'

            if author_on_exec:
                # html += '<img style="width:16px;padding-left:20px;" src="yellow-flag-th.png" alt="Comittee flag" title="At least one author was on the ALSPAC executive
                html += '<div style="text-align:center;font-size:14px;background:#' + config.project_details['colour_hex_secondary'] + ';color:#' + config.project_details['colour_hex_primary'] + ';padding:2px 4px;box-shadow: 0px 0px 1px #4e4e4e inset;">At least one author was a member of the ' + config.project_details['name'] + ' Executive Committee.</div>'
                this_paper['Extras']['author_on_exec'] = True
            else:
                this_paper['Extras']['author_on_exec'] = False

            # Journal volume and issue
            try:
                html += this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
            except:
                pass

            try:
                html += ', Volume ' + this_paper['volume']
            except:
                pass

            try:
                html += ', Issue ' + this_paper['MedlineCitation']['Article']['Journal']['JournalIssue']['Issue']
            except:
                pass
            html += '<br/>'

            # PMID
            try:
                if this_paper['IDs']['PMID']:
                    html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + str(this_paper['IDs']['PMID'])+'">' + str(this_paper['IDs']['PMID']) + '</a>'
            except:
                pass

            # DOI
            try:
                if this_paper['IDs']['DOI']:
                    html += '&nbsp;DOI: <a href="http://doi.org/' + this_paper['IDs']['DOI'] + '">' + this_paper['IDs']['DOI'] + '</a>'
            except:
                pass

            # citation count
            html += "<table class='citation_table'>"
            html += '<tr><th colspan="4">Citation Counts</th></tr>'
            html += '<tr>'
            try:
                html += '<td>Scopus: ' + this_paper['Extras']['Citations'] + '<td>'
            except:
                pass

            try:
                html += '<td>Europe PMC: ' + '-' + '<td>'
            except:
                pass

            html += '</tr>'
            html += "</table>"

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

    # Output the info into an HTML file
    # For each year dict item
    # for this_year in sorted(yearly_papers, key=yearly_papers.get, reverse=True):
    exec_data = {}
    exec_data_string = "year,papers with author on executive committee,total papers"
    for this_year in sorted(yearly_papers, reverse=True):
        # Check there is some data for this year - not all do
        if len(yearly_papers[this_year]) == 0:
            continue

        main += '<a href="' + this_year + '/index.html">' + this_year + '</a> '

        if not os.path.exists(config.html_dir + '/papers/' + this_year):
            os.mkdir(config.html_dir + '/papers/' + this_year)
        year_file = open(config.html_dir + '/papers/' + this_year + '/index.html', 'w')

        shutil.copyfile(config.template_dir + '/altmetric.png', config.html_dir + '/papers/' + this_year + '/altmetric.png')
        shutil.copyfile(config.template_dir + '/yellow-flag-th.png', config.html_dir + '/papers/' + this_year + '/yellow-flag-th.png')

        # Put html together for this page
        temp = '<!DOCTYPE html><html lang="en-GB">'

        # html head
        temp += '<head>'
        temp += '<title>' + site_second_title + '</title>'
        temp += '<link rel="stylesheet" href="../../css/uobcms_corporate.css">'
        temp += '<link rel="stylesheet" href="../../css/style_main.css">'
        temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'
        temp += '</head>'

        temp += build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; <a href="../index.html">Papers List</a> &gt; ' + this_year + '</p>', "../../", "")

        temp += '<h1 id="pagetitle">Papers List - ' + this_year + '</h1>'
        temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

        temp += '<h2>' + str(len(yearly_papers[this_year])) + ' Publications From ' + this_year + '</h2>'
        print >>year_file, temp
        # This is a list
        for this_item in yearly_papers[this_year]:
            temp = this_item.values()
            print >>year_file, temp[0].encode('utf-8')
            for this_paper in papers:
                if this_paper['IDs']['hash'] == this_item.keys()[0]:
                    if this_paper['Extras']['author_on_exec']:
                        try:
                            exec_data[this_year] += 1
                        except:
                            exec_data[this_year] = 1

                    break

            try:
                exec_data[this_year]
            except:
                exec_data[this_year] = 0

        temp = build_common_foot()
        print >>year_file, temp
        try:
            exec_data_string += "\n" + str(this_year) + "," + str(exec_data[this_year]) + "," + str(len(yearly_papers[this_year]))
        except:
            pass

    main += "</p>" + build_common_foot()
    print >>html_file, main

    exec_csv = open(config.data_dir + '/exec_publications.csv', 'w')
    print >>exec_csv, exec_data_string


############################################################
# Build a list of all mesh keywords
############################################################
def build_mesh(papers):

    print "\n###HTML - mesh###"

    shutil.copyfile(config.template_dir + '/keyword_history.js', config.html_dir + '/mesh/keyword_history.js')

    mesh_papers_all = {}
    mesh_papers_major = {}
    html_file_all = open(config.html_dir + '/all_keywords/index.html', 'w')
    html_file_major = open(config.html_dir + '/major_keywords/index.html', 'w')

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
    data_from_count = 0
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

            data_from_count += 1
        except:
            pass

    # Read in mesh tree hierarchy
    f = open(config.config_dir + "/mesh_tree_hierarchy.csv", 'rt')
    mesh_tree = {}
    mesh_tree_reverse = {}
    try:
        reader = csv.reader(f)
        for row in reader:
            mesh_tree[row[2]] = row[0]
            mesh_tree_reverse[row[0]] = row[2]

    finally:
        f.close()

    f = open(config.config_dir + "/mesh_categories.csv", 'rt')
    mesh_categories = {}
    try:
        reader = csv.reader(f)
        for row in reader:
            mesh_categories[row[0]] = row[1]

    finally:
        f.close()

    # Make list of second level mesh headings
    second_found = 0
    top_found = 0
    total = 0
    mesh_second_level_headings = {}
    mesh_top_level_headings = {}
    for this_mesh in mesh_papers_major:
        try:
            tree_number = mesh_tree[this_mesh]
            # Get Second Level
            tree_number_split = tree_number.split(".")
            second_level = mesh_tree_reverse[tree_number_split[0]]

            try:
                mesh_second_level_headings[second_level] += len(mesh_papers_major[this_mesh])
            except:
                mesh_second_level_headings[second_level] = len(mesh_papers_major[this_mesh])
            second_found += len(mesh_papers_major[this_mesh])

            # Get Top Level
            top_level = tree_number[0]
            try:
                mesh_top_level_headings[mesh_categories[top_level]] += len(mesh_papers_major[this_mesh])
            except:
                mesh_top_level_headings[mesh_categories[top_level]] = len(mesh_papers_major[this_mesh])
            top_found += len(mesh_papers_major[this_mesh])
        except:
            pass

        total += len(mesh_papers_major[this_mesh])

    print "MeSH Second Level Found: " + str(second_found) + "/" + str(total)
    print "Unique MeSH Second Level: " + str(len(mesh_second_level_headings))
    print "MeSH Top Level Found: " + str(top_found) + "/" + str(total)
    print "Unique MeSH Top Level: " + str(len(mesh_top_level_headings))

    print "\n" + str(mesh_top_level_headings)
    print "\n" + str(mesh_second_level_headings)

    print "Second Level MeSH"
    for mesh in mesh_second_level_headings:
        print mesh + "\t" + str(mesh_second_level_headings[mesh])

    print "\n\n"

    print "Top Level MeSH"
    for mesh in mesh_top_level_headings:
        print mesh + "\t" + str(mesh_top_level_headings[mesh])

    # Print mesh_papers
    # Make a JSON file for each mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers_all:
        file_name = config.html_dir + '/mesh/all_' + this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers_all[this_mesh], indent=4))
        fo.close()

    # Make a JSON file for each major mesh term, in it put all the PMIDs for this term
    for this_mesh in mesh_papers_major:
        file_name = config.html_dir + '/mesh/major_' + this_mesh
        fo = open(file_name, 'wb')
        fo.write(json.dumps(mesh_papers_major[this_mesh], indent=4))
        fo.close()

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; All Keywords</p>', "../", "")

    temp += '<h1 id="pagetitle">All Keywords</h1>'
    temp += '<p>' + str(len(mesh_papers_all)) + ' Keywords</p>'

    print >>html_file_all, temp

    # Make a page with ALL the headings on it
    print >>html_file_all, '<ul>'
    for this_mesh in sorted(mesh_papers_all):
        temp = '<li><a href="../mesh/' + this_mesh.replace(" ", "%20") + '/index.html">' + this_mesh + '</a></li>'
        print >>html_file_all, temp
    print >>html_file_all, '</ul>'

    temp = build_common_foot()
    print >>html_file_all, temp

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Major Keywords (MeSH)</p>', "../", "")

    temp += '<h1 id="pagetitle">Major Keywords (MeSH)</h1>'
    temp += '<p>' + str(len(mesh_papers_major)) + ' Keywords</p>'

    print >>html_file_major, temp

    # Make a page with the MAJOR headings on it
    print >>html_file_major, '<ul>'
    for this_mesh in sorted(mesh_papers_major):
        temp = '<li><a href="../mesh/' + this_mesh.replace(" ", "%20") + '/index.html">' + this_mesh + '</a></li>'
        print >>html_file_major, temp
    print >>html_file_major, '</ul>'

    temp = build_common_foot()
    print >>html_file_major, temp

    word_cloud_list = "["
    word_cloud_n = 0
    word_cloud_max = 0

    word_cloud_raw = ""

    # MAJOR KEYWORD CLOUD
    # Get first x MAJOR keywords in order
    # Prepare variables
    ordered_mesh_papers_all = {}
    mesh_papers_all_temp = {}
    for this_mesh in mesh_papers_major:
        mesh_papers_all_temp[this_mesh] = mesh_papers_all[this_mesh]

    # Sort and get x words
    for x in range(1, 150):
        max_val = -1
        max_mesh = None
        for this_mesh in mesh_papers_all_temp:
            if len(mesh_papers_all_temp[this_mesh]) > max_val:
                max_val = len(mesh_papers_all_temp[this_mesh])
                max_mesh = this_mesh

        ordered_mesh_papers_all[max_mesh] = mesh_papers_all_temp[max_mesh]
        mesh_papers_all_temp.pop(max_mesh, None)

    # Make papers list for headings pages
    for this_mesh in ordered_mesh_papers_all:

        # Word cloud
        if word_cloud_n > 0:
            word_cloud_list += ','
        word_cloud_n += 1

        number = len(mesh_papers_all[this_mesh])

        if number > word_cloud_max:
            word_cloud_max = number

        word_cloud_list += '{"text":"' + this_mesh + '", "size":' + str(math.sqrt(number)*2.5) + '}'

        for x in range(0, number):
            word_cloud_raw += " " + this_mesh

    # Read in exec csv
    exec_list = []
    f = open(config.config_dir + "/" + config.project_details['short_name'] + '_exec_members.csv', 'rt')
    try:
        reader = csv.reader(f)
        n = 0
        for row in reader:
            if n > 0:
                exec_list.append(row)
            n += 1
    finally:
        f.close()

    # Print page
    for this_mesh in mesh_papers_all:

        if not os.path.exists(config.html_dir + '/mesh/' + this_mesh):
            os.mkdir(config.html_dir + '/mesh/' + this_mesh)

        file_name = config.html_dir + '/mesh/' + this_mesh + '/index.html'
        with codecs.open(file_name, 'wb', "utf-8") as fo:

            # Put html together for this page
            temp = '<!DOCTYPE html><html lang="en-GB">'

            # html head
            temp += '<head>'
            temp += '<title>' + site_second_title + '</title>'
            temp += '<link rel="stylesheet" href="../../css/uobcms_corporate.css">'
            temp += '<link rel="stylesheet" href="../../css/style_main.css">'
            temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'

            temp += '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
            temp += '<script type="text/javascript" src="../' + this_mesh.replace(" ", "%20") + '.js"></script>'
            temp += '<script>var primary_colour = "#' + config.project_details['colour_hex_primary'] + '";</script>'
            temp += '<script>var secondary_colour = "#' + config.project_details['colour_hex_secondary'] + '";</script>'
            temp += '<script type="text/javascript" src="../keyword_history.js"></script>'

            temp += '</head>'

            temp += build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; Keyword &gt; ' + this_mesh + '</p>', "../../", "")

            temp += '<h1 id="pagetitle">Keyword - ' + this_mesh + '</h1>'
            temp += '<h2>All Keyword History</h2>'

            # ===== KEYWORD OVER TIME CALCULATIONS =====
            # First some prep has to be done to set up the array for the number of year. This is copied from the citations graph prep and is probably very inefficent for this task

            summary = {}
            # Calculate the number of papers for each year
            for this_paper in mesh_papers_all[this_mesh]:

                # Get paper object from the hash
                paper_obj = None
                for p in papers:
                    if this_paper == p['IDs']['hash']:
                        paper_obj = p
                        break

                if paper_obj is not None:
                    this_paper = paper_obj
                    try:
                        this_year = this_paper['PubmedData']['History'][0]['Year']
                        # Make sure there is a dict item for this year
                        if this_year not in summary:
                            summary[this_year] = {'num_papers': 0, 'citations': 0}

                        # increment the number of citaitons by one
                        summary[this_year]['num_papers'] += 1

                        # add the citations for this paper to the year running total
                        try:
                            summary[this_year]['citations'] += int(this_paper['Extras']['Citations'])
                        except:
                            pass

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
                        summary[str(this_year)] = {'num_papers': 0, 'citations': 0}

            # Print data to file
            data_file = open(config.html_dir + '/mesh/' + this_mesh + '.js', 'w')
            print >>data_file, 'var papers =([[\'Year\', \'Number of papers\'],'
            for this_year in sorted(summary, reverse=False):
                print >>data_file, '[\''+this_year+'\','+str(summary[this_year]['num_papers'])+'],'
            print >>data_file, ']);'

            print >>data_file, 'var citations =([[\'Year\', \'Number of Citations\'],'
            for this_year in sorted(summary, reverse=False):
                print >>data_file, '[\''+this_year+'\','+str(summary[this_year]['citations'])+'],'
            print >>data_file, ']);'

            temp += '<div id="papers_chart_div"></div>'
            temp += '<div id="citations_chart_div"></div>'

            # List publications
            temp += '<h2>Publications</h2>'
            temp += '<p><em>' + str(len(mesh_papers_all[this_mesh])) + ' publications with this keyword</em></p>'

            print >>fo, temp

            # Build the text needed for each paper
            for this_paper in mesh_papers_all[this_mesh]:

                try:
                    # Get paper object
                    paper_obj = None
                    for p in papers:
                        if this_paper == p['IDs']['hash']:
                            paper_obj = p
                            break

                    if paper_obj is not None:
                        this_paper = paper_obj

                    html = ''

                    # altmetric data
                    try:
                        if this_paper['IDs']['DOI']:
                            html += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>'
                    except:
                        pass

                    # Paper title as a link
                    html += '<span style="text-decoration: underline; font-weight:bold;">' + this_paper['title'] + '</span><br/>'

                    # Authors
                    authors = []
                    author_on_exec = False
                    for this_author in this_paper['author']:
                        # Some author lists have a collective name. Ignore this.
                        # Some people don't actually have initials. eg wraight in pmid:18454148
                        try:
                            # Check if an author was on the exec Comittee
                            for x in exec_list:
                                # Check if authors name matches
                                if x[2] == this_author['clean']:
                                    # Get start and end date of exec membership
                                    exec_start = time.mktime(datetime.datetime.strptime(x[0], "%d/%m/%Y").timetuple())
                                    exec_end = int(time.time())
                                    if not x[1] == "":
                                        exec_end = time.mktime(datetime.datetime.strptime(x[1], "%d/%m/%Y").timetuple())

                                    # Convert issued date into a timestamp
                                    issued_date = this_paper['issued']['date-parts']
                                    date = str(issued_date[0][2]) + "/" + str(issued_date[0][1]) + "/" + str(issued_date[0][0])
                                    issued_timestamp = time.mktime(datetime.datetime.strptime(date, "%d/%m/%Y").timetuple())

                                    # If publication is issued between exec_start and exec_end then flag
                                    if issued_timestamp > exec_start and issued_timestamp < exec_end:
                                        author_on_exec = True
                                        break

                            authors.append(this_author['family'] + ', ' + this_author['given'])
                        except:
                            pass

                    html += htmlentities.encode('; '.join(authors))
                    html += '<br/>'

                    if author_on_exec:
                        # html += '<img style="width:16px;padding-left:20px;" src="yellow-flag-th.png" alt="Comittee flag" title="At least one author was on the ALSPAC executive
                        html += '<div style="text-align:center;font-size:14px;background:#' + config.project_details['colour_hex_secondary'] + ';color:#' + config.project_details['colour_hex_primary'] + ';padding:2px 4px;box-shadow: 0px 0px 1px #4e4e4e inset;">At least one author was a member of the ' + config.project_details['name'] + ' Executive Committee.</div>'
                        this_paper['Extras']['author_on_exec'] = True
                    else:
                        this_paper['Extras']['author_on_exec'] = False

                    # Journal volume and issue
                    try:
                        html += this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
                    except:
                        pass

                    try:
                        html += ', Volume ' + this_paper['volume']
                    except:
                        pass

                    try:
                        html += ', Issue ' + this_paper['MedlineCitation']['Article']['Journal']['JournalIssue']['Issue']
                    except:
                        pass
                    html += '<br/>'

                    # PMID
                    try:
                        if this_paper['IDs']['PMID']:
                            html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + str(this_paper['IDs']['PMID'])+'">' + str(this_paper['IDs']['PMID']) + '</a>'
                    except:
                        pass

                    # DOI
                    try:
                        if this_paper['IDs']['DOI']:
                            html += '&nbsp;DOI: <a href="http://doi.org/' + this_paper['IDs']['DOI'] + '">' + this_paper['IDs']['DOI'] + '</a>'
                    except:
                        pass

                    # citation count
                    html += "<table class='citation_table'>"
                    html += '<tr><th colspan="4">Citation Counts</th></tr>'
                    html += '<tr>'
                    try:
                        html += '<td>Scopus: ' + this_paper['Extras']['Citations'] + '<td>'
                    except:
                        pass

                    try:
                        html += '<td>Europe PMC: ' + '-' + '<td>'
                    except:
                        pass

                    html += '</tr>'
                    html += "</table>"

                    # Add an extra line break at the end
                    html += '<br/><br/>'

                    fo.write(html)

                except:
                    pass

            temp = build_common_foot()
            print >>fo, temp

        fo.close()

    word_cloud_list += "]"
    build_word_cloud(papers, word_cloud_list, data_from_count)


###########################################################
# Build a google map based on the lat longs provided before.
###########################################################
def build_google_map(papers):

    print "\n###HTML - Insititutions Map###"

    info = []
    number_of_points = 0
    for this_paper in papers:
        try:
            this_place = {'lat': this_paper['Extras']['LatLong']['lat'], 'long': this_paper['Extras']['LatLong']['long'], 'name': this_paper['Extras']['CleanInstitute']}
            info.append(this_place)
            number_of_points += 1
        except:
            pass

    # print info

    kml = "var locations =["
    for this_info in info:
        kml += '["' + this_info['name'] + '",' + str(this_info['lat']) + ',' + str(this_info['long']) + '],'
    kml += ']'

    kml_file = codecs.open(config.html_dir + '/map/map.kml', 'w', 'utf-8')
    print >>kml_file, kml

    html_file = codecs.open(config.html_dir + '/map/index.html', 'w', 'utf-8')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="http://maps.googleapis.com/maps/api/js?sensor=false"></script>'
    temp += '<script type="text/javascript" src="map.kml"></script>'
    temp += '<script type="text/javascript" src="map.js"></script>'

    shutil.copyfile(config.template_dir + '/map.js', config.html_dir + '/map/map.js')
    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/map/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/css/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Institutions Map</p>', "../", "onload='initialize()'")

    temp += '<h1 id="pagetitle">Institutions Map</h1>'

    temp += "<div class='loading'><img src='loading.gif' alt='Loading'></div>"
    temp += "<div id='map_canvas'></div>"
    temp += "<p>Data from " + intWithCommas(number_of_points) + " publications</p>"

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# Publications by country
###########################################################
def build_country_map(papers, api_key):

    print "\n###HTML - Country Map###"

    countries = {}
    number_of_points = 0
    for this_paper in papers:
        try:

            if this_paper['Extras']['country_code'] in countries:
                countries[this_paper['Extras']['country_code']] += 1
            else:
                countries[this_paper['Extras']['country_code']] = 1
            number_of_points += 1
        except:
            pass

    country_string = ""
    for country in countries.keys():
        country_string += ",['" + country + "'," + str(countries[country]) + "]"

    html_file = open(config.html_dir + '/country/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="https://maps.googleapis.com/maps/api/js?key=' + config.google_maps_api_key + '&libraries=visualization"></script>'
    temp += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script> <script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    temp += '<script type="text/javascript">' + "google.charts.load('current', {'packages':['geochart']});google.charts.setOnLoadCallback(drawRegionsMap);function drawRegionsMap() {var data = google.visualization.arrayToDataTable([ ['Country', 'Publications']" + country_string + "]); var options = { colorAxis: {colors: ['#" + config.project_details['colour_hex_secondary'] + "', '#" + config.project_details['colour_hex_primary'] + "']} }; var chart = new google.visualization.GeoChart(document.getElementById('regions_div')); chart.draw(data, options); }</script>"

    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/country/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/country/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications by Country</p>', "../", "")

    temp += '<h1 id="pagetitle">Publications by Country</h1>'

    temp += "<div class='loading'><img src='loading.gif' alt='Loading'></div>"
    temp += "<div id='regions_div' style='width: 900px; height: 500px;'></div>"
    temp += "<p>Data from " + intWithCommas(number_of_points) + " publications</p>"

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp

    build_city_map(papers)


###########################################################
# Publications by UK city
###########################################################
def build_city_map(papers):

    print "\n###HTML - City Map###"

    cities = {}
    number_of_points = 0
    for this_paper in papers:

        try:
            if this_paper['Extras']['postal_town'] != "":
                if this_paper['Extras']['postal_town'] in cities:
                    cities[this_paper['Extras']['postal_town']] += 1
                else:
                    cities[this_paper['Extras']['postal_town']] = 1
                number_of_points += 1
        except:
            pass

    city_string = ""
    for city in cities.keys():
        city_string += ",['" + city + "'," + str(cities[city]) + "]"

    html_file = open(config.html_dir + '/city/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script> <script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    temp += "<script>google.charts.load('current', {'packages':['geochart']});google.charts.setOnLoadCallback(drawMarkersMap);function drawMarkersMap() {var data = google.visualization.arrayToDataTable([['City',   'Publications']" + city_string + " ]); var options = {region: 'GB', displayMode: 'markers', colorAxis: {colors: ['#" + config.project_details['colour_hex_secondary'] + "', '#" + config.project_details['colour_hex_primary'] + "']}}; var chart = new google.visualization.GeoChart(document.getElementById('regions_div'));chart.draw(data, options); };</script>"

    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/city/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/city/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications by UK City</p>', "../", "")

    temp += '<h1 id="pagetitle">Publications by UK City</h1>'

    temp += "<div class='loading'><img src='loading.gif' alt='Loading'></div>"
    temp += '<div id="regions_div" style="width: 900px; height: 500px;"></div>'
    temp += "<p>Data from " + intWithCommas(number_of_points) + " publications</p>"
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


def build_metrics(papers, cohort_rating, cohort_rating_data_from, study_start_year, study_current_year):

    print "\n###HTML - Metrics###"

    html_file = open(config.html_dir + '/metrics/index.html', 'w')

    # NUMBER OF PAPERS PER CITATION COUNT
    num_papers_citations = []
    max_citations = 0
    for this_paper in papers:
        try:
            n_cits = int(this_paper['Extras']['Citations'])
            if n_cits > max_citations and n_cits < 200:
                max_citations = n_cits
            try:
                num_papers_citations[n_cits] += 1
            except:
                num_papers_citations.insert(n_cits, 1)
        except:
            pass

    # Create data string for plot
    n_papers_with_x_citations = "var papers_per_citation_count = ([['Number of Citations','Number of Papers']"
    # Add Zeros in missing indexes
    for this_n_citations in range(1, max_citations+1):
        try:
            n_papers_with_x_citations += ",[" + str(this_n_citations) + "," + str(num_papers_citations[this_n_citations]) + "]"
        except:
            n_papers_with_x_citations += ",[" + str(this_n_citations) + ",0]"

    n_papers_with_x_citations += "])"

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'

    shutil.copyfile(config.template_dir + '/metrics.js', config.html_dir + '/metrics/metrics.js')

    temp += '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    temp += '<script type="text/javascript" src="../data.js"></script>'
    temp += '<script>' + n_papers_with_x_citations + '</script>'
    temp += '<script type="text/javascript" src="../map/map.js"></script>'
    temp += '<script>var primary_colour = "#' + config.project_details['colour_hex_primary'] + '";</script>'
    temp += '<script type="text/javascript" src="metrics.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Study Metrics</p>', "../", "")

    temp += '<h1 id="pagetitle">Study Metrics</h1>'

    temp += "<p>{Explanation of study metrics}</p>"

    # Metric calculations
    total_publications = len(papers)
    total_citations = 0
    total_citations_data_from_count = 0
    paper_citations = []
    c20_index = 0

    # study_start_year = 1991
    # study_current_year = 2016
    study_duration = study_current_year - study_start_year

    c_index_bound = 100

    for this_paper in papers:
        # add the citations for paper to total
        try:

            cit = int(this_paper['Extras']['Citations'])
            total_citations += cit
            total_citations_data_from_count += 1
            paper_citations.append(cit)

            # increment c20-index if more that 20 citations
            if cit >= c_index_bound:
                c20_index += 1
        except:
            pass

    average_citations = float(total_citations)/float(total_citations_data_from_count)
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
    temp += "<div class='metric_name'>h-index</div>"
    temp += "<div class='metric_value'>" + str(h_index) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>h-index is the largest number h such that h publications from a study have at least h citations.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>g-index</div>"
    temp += "<div class='metric_value'>" + str(g_index) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The largest number n of highly cited articles for which the average number of citations is at least n.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Cohort-Rating</div>"
    temp += "<div class='metric_value'>" + str("{0:.3f}".format(round(cohort_rating, 3))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(cohort_rating_data_from) + " Publications</div>"
    temp += "<div class='metric_description'></div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Mean Citations Per Publication</div>"
    temp += "<div class='metric_value'>" + str("{0:.2f}".format(round(average_citations, 2))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The total number of citations divided by the total number of publications.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>c" + str(c_index_bound) + "-index</div>"
    temp += "<div class='metric_value'>" + intWithCommas(c20_index) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The number of publications from a study that have at least " + str(c_index_bound) + " citations.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>c" + str(c_index_bound) + "-index per Study Year</div>"
    temp += "<div class='metric_value'>" + str("{0:.2f}".format(round(i20_index_per_year, 2))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The c" + str(c_index_bound) + "-index divided by the number of years the study has been running for.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Total Publications</div>"
    temp += "<div class='metric_value'>" + intWithCommas(total_publications) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_publications) + " Publications</div>"
    temp += "<div class='metric_description'>This is the number of all publications for the study.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Total Citations</div>"
    temp += "<div class='metric_value'>" + intWithCommas(total_citations) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>This is the number of citations to all publications for the study.</div>"
    temp += "</div>"
    temp += "</div>"
    temp += "<div class='clear'></div>"

    temp += '<div id="cumulative_div"></div>'
    temp += "<p style='text-align:center;'>Data from " + intWithCommas(cohort_rating_data_from) + " publications</p>"
    temp += '<div id="papers_per_year_div"></div>'
    temp += "<p style='text-align:center;'>Data from " + intWithCommas(cohort_rating_data_from) + " publications</p>"
    temp += '<div id="papers_per_citation_count_div"></div>'
    temp += "<p style='text-align:center;'>Data from " + intWithCommas(cohort_rating_data_from) + " publications</p>"
    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# Build keyword word cloud
###########################################################
def build_word_cloud(papers, list, data_from_count):

    print "\n###HTML - Keyword Cloud###"

    html_file = open(config.html_dir + '/wordcloud/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'
    temp += '<style>.wordcloud{ width:100%; height:500px;}</style>'

    shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/wordcloud/d3wordcloud.js')
    shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/wordcloud/d3.layout.cloud.js')

    temp += '<script>var word_list = ' + list + '</script>'
    temp += '<script src="http://d3js.org/d3.v3.min.js"></script>'
    temp += '<script src="d3.layout.cloud.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Major Keyword Cloud</p>', "../", "")

    temp += '<h1 id="pagetitle">Major Keyword Cloud</h1>'

    temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
    temp += '</cloud>'
    temp += "<p>Data from " + intWithCommas(data_from_count) + " publications</p>"

    temp += '<script src="d3wordcloud.js"></script>'

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# Build abstract word cloud
###########################################################
def build_abstract_word_cloud(papers, data_from_count):

    print "\n###HTML - Abstract Word Cloud###"

    f = open(config.data_dir + "/abstracts.csv", 'rt')

    list = "["
    n = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            if n > 0:
                list += ","

            if row[0] != "":
                # list += '["' + row[0].replace("'","\'").replace('"','\"') + '",' + str(row[1]) +  ']'
                list += '{"text":"' + row[0].replace("'", "\'").replace('"', '\"') + '","size":' + str(math.sqrt(int(row[1]))*1.5) + '}'
                n += 1

            if n > 200:
                break

    finally:
        f.close()

    list += "];"
    list_file = open(config.html_dir + '/abstractwordcloud/list.js', 'w')
    print >>list_file, " var word_list = " + list

    html_file = open(config.html_dir + '/abstractwordcloud/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'
    temp += '<style>.wordcloud{ width:100%; height:500px;}</style>'

    shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/abstractwordcloud/d3wordcloud.js')
    shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/abstractwordcloud/d3.layout.cloud.js')

    temp += '<script src="list.js"></script>'
    temp += '<script src="http://d3js.org/d3.v3.min.js"></script>'
    temp += '<script src="d3.layout.cloud.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Abstract Word Cloud</p>', "../", "")

    temp += '<h1 id="pagetitle">Abstract Word Cloud</h1>'

    temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
    temp += '</cloud>'

    temp += '<script src="d3wordcloud.js"></script>'
    temp += "<p>Data from " + intWithCommas(data_from_count) + " publications</p>"

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# Build Author Network
###########################################################
def get_author_string_from_hash(hash_string, network):

    for author in network['authors']:
        if author == hash_string:
            return network['authors'][author]['clean']


def build_author_network(papers, network):

    print "\n###HTML - Author Network###"

    # Create json file
    net_file = open(config.html_dir + '/authornetwork/network.json', 'w')

    net_json = '{'
    net_json += '"nodes":['
    n = 0
    print >>net_file, net_json

    for author in network['authors']:
        # print network['authors'][author]
        net_json = ""
        if n > 0:
            net_json += ","
        net_json += '{"id": "' + network['authors'][author]['clean'] + '", "group":1}'

        print >>net_file, net_json
        n += 1

    net_json = '],"links": ['
    print >>net_file, net_json

    n = 0
    for con in network['connections']:
        try:
            net_json = ""
            if n > 0:
                net_json += ","

            author_0 = get_author_string_from_hash(network['connections'][con]['authors'][0]['author_hash'], network)
            author_1 = get_author_string_from_hash(network['connections'][con]['authors'][1]['author_hash'], network)

            n_con = network['connections'][con]['num_connections']/2

            net_json += '{"source": "' + author_0 + '", "target": "' + author_1 + '", "value": ' + str(n_con) + '}'

            print >>net_file, net_json
            n += 1
        except:
            pass

    net_json = ']'
    net_json += '}'
    print >>net_file, net_json

    html_file = open(config.html_dir + '/authornetwork/index.html', 'w')

    shutil.copyfile(config.template_dir + '/network.js', config.html_dir + '/authornetwork/network.js')
    shutil.copyfile(config.template_dir + '/author_network.png', config.html_dir + '/authornetwork/author_network.png')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<style>.links line {  stroke: #999;  stroke-opacity: 0.6;} .nodes circle {  stroke: #fff;  stroke-width: 1.5px;} </style>'

    temp += '<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.1.0/jquery.min.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Author Network</p>', "../", "")

    temp += '<h1 id="pagetitle">Author Network</h1>'

    # Print nodes to csv
    nodes_csv = open(config.html_dir + '/authornetwork/nodes.csv', 'w')

    print >>nodes_csv, 'id,Label'
    n = 0

    for author in network['authors']:
        print >>nodes_csv, author + "," + network['authors'][author]['clean']
        n += 1

    # Print conections to csv
    connections_csv = open(config.html_dir + '/authornetwork/connections.csv', 'w')

    print >>connections_csv, 'Source,Target'

    n = 0
    for con in network['connections']:
        try:

            author_0 = network['connections'][con]['authors'][0]['author_hash']
            author_1 = network['connections'][con]['authors'][1]['author_hash']

            n_con = network['connections'][con]['num_connections']/2

            print >>connections_csv, '"' + author_0 + '","' + author_1 + '"'

        except:
            pass
        n += 1

    temp += '<a href="author_network.png"><img src="author_network.png" alt="author network"></a>'

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# ErrorLog
###########################################################
def build_error_log(papers, error_log):

    print "\n###HTML - Error Log###"

    html_file = open(config.html_dir + '/errorlog/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script>function jsonFormat(){ var x = document.getElementsByClassName("textoutput"); for (i = 0; i < x.length; i++) { x[i].innerHTML = JSON.stringify(JSON.parse(x[i].innerHTML) , null, 4); } }</script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Error Log</p>', "../", "onload='jsonFormat();'")

    temp += '<h1 id="pagetitle">Error Log</h1>'

    temp += error_log.printLog()

    temp += "<script>"
    # Set cookie
    temp += ' var expiration_date = new Date(); var cookie_string = ""; expiration_date.setFullYear(expiration_date.getFullYear() + 1); cookie_string = "show_error_page=true; path=/; expires=" + expiration_date.toUTCString(); document.cookie = cookie_string;'
    temp += "</script>"

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# Help Page
###########################################################
def build_help():

    print "\n###HTML - Help/Information page###"

    html_file = open(config.html_dir + '/help/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/uobcms_corporate.css">'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Information</p>', "../", "")

    temp += '<h1 id="pagetitle">Information</h1>'

    temp += '<h2>Where does the data come from?</h2>'
    temp += '<h3>Publication Data</h3>'
    temp += '<p>Data and metadata for the publications are located using the PubMed search engine.</p>'

    temp += '<h3>Citations Data</h3>'
    temp += '<p>Citation data is retrieved from the Scopus API provided by Elsevier.</p>'

    temp += '<h3>Geodata</h3>'
    temp += '<p>The Coordinate location of institutions is retrieved from the free project Wikidata. The coordinate data is then converted into country'
    temp += ' and city names using the Google Maps API.</p>'

    temp += "<h2>Why don't some statistics use data from all publications?</h2>"

    temp += '<p>Data can be missing because some publications are very old. There are many different ways to track publications'
    temp += ' and Prior to the introduction of DOIs in 2000 there was no unified method. This means some metadata on old publications could have been lost or not recorded.</p>'

    temp += '<p>The data used for the statistics are gathered from databases which only collect data from particular journals.'
    temp += ' This problem is most obvious for citation data from Scopus. Although they have a particular publication on record,'
    temp += ' there may be citations for this publication from a journal that they do not index and therefore these will not be in the citation count.</p>'

    print >>html_file, temp

    temp = build_common_foot()
    print >>html_file, temp


###########################################################
# CSS colour scheme
###########################################################
def build_css_colour_scheme():

    html_file = open(config.template_dir + '/colour_scheme.css', 'w')

    temp = ".uob-header-container {background: #" + config.project_details['colour_hex_primary'] + ";}"

    temp += "/* h1 */"
    temp += "#pagetitle {"
    temp += "color:#" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    temp += ".maintitle a {"
    temp += "color: #fff!important;"
    temp += "}"

    temp += ".maintitle {"
    temp += "border-left: 1px solid rgba(255, 255, 255, 0.2);"
    temp += "}"

    temp += "#uoblogo span {"
    temp += "background-image: url(" + config.project_details['header_image_url'] + ");"
    temp += "}"

    temp += ".uob-header-container:before,"
    temp += ".uob-header-container:after {"
    temp += "display: block;"
    temp += "visibility: visible;"
    temp += "}"

    temp += ".nav-users a {"
    temp += "color: #fff;"
    temp += "}"

    temp += "color: rgba(255,255,255,0.4)!important;"
    temp += ".nav-users .icon-arrow-right:after {"
    temp += "}"

    temp += "/* keywords list */"
    temp += "div#uobcms-content ul a {"
    temp += "color: #" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    temp += "/* Metrics Page */"

    temp += "div#uobcms-content div.metric {"
    temp += "background:#efede9;"
    temp += "color:#" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    temp += "div#uobcms-content div.metric:hover {"
    temp += "background:#" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    temp += "div#uobcms-content div.metric:hover div.metric_name {"
    temp += "color:#f2f2f2;"
    temp += "}"

    temp += "div#uobcms-content div.metric:hover div.metric_value {"
    temp += "color:#f2f2f2;"
    temp += "}"

    temp += "div#uobcms-content div.metric div.metric_description {"
    temp += "color:#efede9;"
    temp += "}"

    temp += "a:link {color:#" + config.project_details['colour_hex_primary'] + "}"
    temp += "a:visited {color:#" + config.project_details['colour_hex_primary'] + "}"

    print >>html_file, temp
