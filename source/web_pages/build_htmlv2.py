#!/usr/bin/env python3

import json
import shutil
import os.path
import csv
from html import escape
import time
import codecs
import os
import logging
import sys
import math

import config.config as config

# The colour scheme is automatic. The colours from the config are automatically put into the pages.
# This .py file also handles the generation of some CSS files.

site_second_title = " Data Set Publications"


############################################################
# Common Page Features
############################################################
def build_common_body(breadcrumb, nav_path, body):
    # This function builds the common header and nav bar for all pages.
    # nav_path used for changes to relative pathing depending on the page (i.e. Home does not need anything but next level down needs leading ../)
    # body is used for putting extra attributes into the body tag (e.g. onload="function();")

    html = "<body " + body + ">"

    html += "<div class='header-container'>"
    html += "<div class='header width-master' role='banner'>"

    if os.path.isfile(config.config_dir + '/' + config.project_details['header_institution_logo_filename']):
        shutil.copy(config.config_dir + '/' + config.project_details['header_institution_logo_filename'], config.html_dir + '/' + config.project_details['header_institution_logo_filename'])
        html += '<div id="main-logo"><a accesskey="1" title="' + config.project_details['header_institution_name'] + '" href="' + config.project_details['header_institution_url'] + '"><img src="' + nav_path + config.project_details['header_institution_logo_filename'] + '" alt="Institution logo"/></a></div>'

    html += "<div class='maintitle' id='maintitle1'>"
    html += "<span id='title1'><a href='" + nav_path + "index.html'>" + config.project_details['name'] + " - " + site_second_title + "</a></span>"
    html += "</div>"
    html += "</div>"
    html += "</div>"

    html += '<div id="wrapper" class="width-master">'
    html += '<div id="col1" role="navigation">'
    html += '<!-- navigation object : navigation -->'

    html += '<ul class="navgroup">'
    html += '<li><a href="' + nav_path + 'index.html">Home</a></li>'
    html += '<li><a href="' + nav_path + 'help/index.html">About</a></li>'
    html += '<li><a href="' + nav_path + 'search/index.html">Search</a></li>'
    html += '<li><a href="' + nav_path + 'all_keywords/index.html">All Keywords</a></li>'
    html += '<li><a href="' + nav_path + 'major_keywords/index.html">Major Keywords (MeSH)</a></li>'
    html += '<li><a href="' + nav_path + 'country/index.html">Map by Country</a></li>'
    html += '<li><a href="' + nav_path + 'institute/index.html">Map by UK institute</a></li>'

    if config.page_show_author_network == "True":
        html += '<li><a href="' + nav_path + 'authornetwork/index.html">Author Network</a></li>'

    html += '<li><a href="' + nav_path + 'metrics/index.html">Study Metrics</a></li>'
    html += '<li><a href="' + nav_path + 'keyword_wordcloud/index.html">Keyword Cloud</a></li>'
    html += '<li><a href="' + nav_path + 'abstractwordcloud/index.html">Abstract Word Cloud</a></li>'

    if not config.public_facing:
        html += '<li><a href="' + nav_path + 'coverage_report.html">Coverage report</a></li>'

    html += '</ul>'

    html += '<div class="after-navgroup">'
    html += '<!-- navigation object : navigation bottom -->'
    html += '<!-- start navigation : additional logo -->'

    # Add the side logo to the actual project if it has been set
    if os.path.isfile(config.config_dir + '/' + config.project_details['side_image_filename']):
        shutil.copy(config.config_dir + '/' + config.project_details['side_image_filename'], config.html_dir + '/' + config.project_details['side_image_filename'])
        html += '<div class="logo-additional">'
        html += '<a href="' + config.project_details['side_image_link'] + '"><img style="max-height:200px" src="' + nav_path + config.project_details['side_image_filename'] + '" alt=""/></a>'
        html += '</div>'

    html += '</div>'
    html += '</div>'

    html += breadcrumb

    html += '<div id="content" role="main">'

    return html


###########################################################
# Build a common footer for all the pages.
###########################################################
def build_common_foot():

    # Output the timestamp to file
    timestamp_file = open(config.html_dir + '/timestamp.txt', 'w')
    timestamp_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))
    timestamp_file.close()

    html = '</div>'
    html += '</div>'
    html += '<div class="foot">'
    # Put a link to the timestamp file, but update it to the content via JS. This means
    # there will be a mechanism to see this when running locally as COS stops JS calls.
    html += '<span id="update_timestamp"><a href="' + config.html_dir + '/timestamp.html">Update time</a></span>'
    html += "<script type='text/javascript' src='../timestamp.js'></script>"
    html += '<script>update_timestamp("/puma/' + config.project_details['short_name'] + '/timestamp.txt")</script>'
    # html += ' Stats last updated on ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
    html += '</div>'
    html += '</body>'
    html += '</html>'

    return html


############################################################
# Home page with summary of years
############################################################
def build_home(papers):

    print("\n###HTML - Home###")

    summary = {}
    missing_year = {'num_papers': 0, 'citations': 0}

    html_file = open(config.html_dir + '/index.html', 'w')
    data_file = open(config.html_dir + '/data.js', 'w')

    # Copy CSS files
    shutil.copyfile(config.template_dir + '/style_main.css', config.html_dir + '/css/style_main.css')

    shutil.copyfile(config.template_dir + '/timestamp.js', config.html_dir + '/timestamp.js')


    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="css/style_main.css">'
    temp += '<link rel="stylesheet" href="css/colour_scheme.css">'
    temp += '<script src="timestamp.js"></script>'
    temp += '</head>'

    temp += build_common_body("", "", "")
    temp += '<h1 id="pagetitle">Summary by Year</h1>'

    html_file.write(temp)

    # Calculate the number of papers for each year
    for this_paper in papers:
        try:
            this_year = this_paper['clean']['clean_date']['year']

            # Make sure there is a dict item for this year
            if this_year not in summary:
                summary[this_year] = {'num_papers': 0, 'cumulative': 0, 'citations': 0, 'cumulative_citations': 0}

            # increment the number of citations by one
            summary[this_year]['num_papers'] += 1

            # add the citations for this paper to the year running total
            try:
                summary[this_year]['citations'] += int(this_paper['clean']['citations']['scopus']['count'])
            except:
                pass

        except:
            try:
                this_paper['clean']['clean_date']['year']
            except:
                missing_year['num_papers'] += 1
                try:
                    missing_year['citations'] += int(this_paper['clean']['citations']['scopus']['count'])
                except:
                    pass

    # Add in some zeros when there is no papers for this year
    years = list(summary.keys())
    first_year = min(years)
    last_year = max(years)
    for this_year in range(int(first_year), int(last_year)):
        try:
            summary[str(this_year)]['num_papers']
        except:
            summary[str(this_year)] = {'num_papers': 0, 'cumulative': 0, 'citations': 0, 'cumulative_citations': 0}

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
    data_file.write('var cumulative =([[\'Year\', \'Number of papers\'],')
    for this_year in sorted(summary, reverse=False):
        data_file.write('[\''+this_year+'\','+str(summary[this_year]['cumulative']) + '],')
    data_file.write(']);')

    # Number per year now
    data_file.write('var papers_per_year=([[\'Year\', \'Number of papers\'],')
    for this_year in sorted(summary, reverse=False):
        data_file.write('[\''+this_year+'\','+str(summary[this_year]['num_papers']) + '],')
    data_file.write(']);')

    # Cohort age weighted citation calculation
    cr_current_year = float(config.metrics_study_current_year)
    cr_sum = 0.0
    cr_data_from = 0

    # print summary
    # Make a page with the headings on it
    html_file.write('<table>')
    html_file.write('<tr><th>Year</th><th>Number published</th><th>Cumulative</th><th>Citations* for papers published in this year</th><th>Cumulative citations* for papers published in this year</th></tr>')
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
        temp += '<td>' + intWithCommas(summary[this_year]['citations']) + '</td>'
        temp += '<td>' + intWithCommas(summary[this_year]['cumulative_citations']) + '</td></tr>'
        html_file.write(temp)

    # Unknown Row
    if missing_year['num_papers'] > 0:
        temp = '<tr>'
        temp += '<td style="font-size:12px;font-weight:bold;"><a href="papers/unknown/index.html">UNKNOWN</a></td>'
        temp += '<td>' + intWithCommas(missing_year['num_papers']) + '</td>'
        temp += '<td>-</td>'
        temp += '<td>' + intWithCommas(missing_year['citations']) + '</td>'
        temp += '<td>-</td>'
        temp += '</tr>'
        html_file.write(temp)
    html_file.write('</table>')

    temp = "<p>Publication year known for " + intWithCommas(cr_data_from) + " of " + intWithCommas(len(papers)) + " publications</p>"
    temp += '<p>* Citation data from <a href="https://www.scopus.com">Scopus</a>.</p>'

    temp += build_common_foot()
    html_file.write(temp)

    cr_sum = cr_sum / len(papers)
    return cr_sum, cr_data_from


############################################################
# Draw Paper Function
# This function is used in the different paper lists to
# display a consistently formatted list.
############################################################
def draw_paper(this_paper):
    html = '<div class="paper">'

    # altmetric data
    try:
        if this_paper['IDs']['DOI']:
            html += '<div style="float:right;" data-badge-popover="right" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>'
    except:
        pass

    # Paper title
    html += '<span style="text-decoration: underline; font-weight:bold;">' + this_paper['clean']['title'] + '</span><br/>'

    # Authors
    authors = []
    for this_author in this_paper['clean']['full_author_list']:
        # Some author lists have a collective name. Ignore this.
        # Some people don't actually have initials. eg wraight in pmid:18454148
        try:
            authors.append(this_author['family'] + ', ' + this_author['given'])
        except:
            logging.debug('This author dropped from author list on webpage for some reason: ' + str(this_author))

    html += escape('; '.join(authors))
    html += '<br/>'

    html += this_paper['clean']['journal']['journal_name']

    try:
        html += ', Volume ' + this_paper['clean']['journal']['volume']
    except:
        pass

    try:
        html += ', Issue ' + this_paper['clean']['journal']['issue']
    except:
        pass
    html += '<br/>'

    # PMID
    try:
        if this_paper['IDs']['PMID']:
            html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + str(this_paper['IDs']['PMID'])+'">' + str(this_paper['IDs']['PMID']) + '</a>&nbsp;'
    except:
        pass

    # DOI
    try:
        if this_paper['IDs']['DOI']:
            html += 'DOI: <a href="https://doi.org/' + this_paper['IDs']['DOI'] + '">' + this_paper['IDs']['DOI'] + '</a>&nbsp;'
    except:
        pass

    # Citation Counts and Sources
    number_citations_counts = 1  # The number of different citation count sources
    citations_counts_width = 100 / number_citations_counts
    html += "<table class='citation_table'>"
    html += '<tr><th colspan="' + str(number_citations_counts) + '">Citation Counts</th></tr>'
    html += '<tr>'
    try:
        # Try to display citation count with link to scopus page
        scopus_links = this_paper['raw']['scopus_data']['search-results']['entry'][0]['link']
        for this_link in scopus_links:
            if this_link['@ref'] == 'scopus-citedby':
                html += '<td style="width:' + str(citations_counts_width) + '%;">Scopus: <a href="' + this_link['@href'] + '" target="_blank">' + str(this_paper['clean']['citations']['scopus']['count']) + '</a></td>'

        # html += '<td style="width:' + str(citations_counts_width) + '%;">Scopus: <a href="https://www.scopus.com/record/display.uri?eid=' + str(this_paper['IDs']['scopus']) + '&origin=inward&txGid=0">' + str(this_paper['clean']['citations']['scopus']['count']) + '</a></td>'
    except:
        try:
            html += '<td style="width:' + str(citations_counts_width) + '%;">Scopus: ' + str(this_paper['clean']['citations']['scopus']['count']) + '</td>'
        except:
            html += '<td style="width:' + str(citations_counts_width) + '%;">Scopus: -</td>'

    html += '</tr>'
    html += "</table>"
    html += '</div>'

    return html


############################################################
# Home page with summary of years
############################################################
def build_papers(papers):

    print("\n###HTML papers list###")

    yearly_papers = {}
    html_file = open(config.html_dir + '/papers/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Papers List</p>', "../", "")

    temp += '<h1 id="pagetitle">Papers List</h1>'
    # Altmetric include
    temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

    html_file.write(temp)
    main = "<p>"

    # Build the text needed for each paper
    for this_paper in papers:
        try:
            # Call draw paper function
            html = draw_paper(this_paper)

            # Append this paper to the list indexed by the year
            this_year = this_paper['clean']['clean_date']['year']

            # Make sure there is a dict item for this year
            if this_year not in yearly_papers:
                yearly_papers[this_year] = list()

            temp = yearly_papers[this_year]
            temp.append({this_paper['IDs']['hash']: html})
            yearly_papers[this_year] = temp
        except:
            print('Failing on ' + this_paper['IDs']['hash'])
            print(sys.exc_info())

    # Output the info into an HTML file
    # For each year dict item
    for this_year in sorted(yearly_papers, reverse=True):
        # Check there is some data for this year - not all do
        if len(yearly_papers[this_year]) == 0:
            continue

        main += '<a href="' + this_year + '/index.html">' + this_year + '</a><br/> '

        if not os.path.exists(config.html_dir + '/papers/' + this_year):
            os.mkdir(config.html_dir + '/papers/' + this_year)
        year_file = open(config.html_dir + '/papers/' + this_year + '/index.html', 'w')

        # Put html together for this page
        temp = '<!DOCTYPE html><html lang="en-GB">'

        # html head
        temp += '<head>'
        temp += '<title>' + site_second_title + '</title>'
        temp += '<link rel="stylesheet" href="../../css/style_main.css">'
        temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'
        temp += '</head>'

        temp += build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; <a href="../index.html">Papers List</a> &gt; ' + this_year + '</p>', "../../", "")

        temp += '<h1 id="pagetitle">Papers List - ' + this_year + '</h1>'
        temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

        temp += '<h2>' + str(len(yearly_papers[this_year])) + ' Publications From ' + this_year + '</h2>'
        year_file.write(temp)
        # This is a list
        for this_item in yearly_papers[this_year]:
            temp = list(this_item.values())
            year_file.write(temp[0])

        temp = build_common_foot()
        year_file.write(temp)

    main += "</p>" + build_common_foot()
    html_file.write(main)

    # == Publications from unknown years ==
    if not os.path.exists(config.html_dir + '/papers/unknown'):
        os.mkdir(config.html_dir + '/papers/unknown')
    unknown_file = open(config.html_dir + '/papers/unknown/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; <a href="../index.html">Papers List</a> &gt; Unknown Year</p>', "../../", "")

    temp += '<h1 id="pagetitle">Papers List - Unknown Year</h1>'
    # Altmetric include
    temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

    n = 0
    html = ""
    for this_paper in papers:
        try:
            this_paper['clean']['clean_date']['year']
        except:
            html += draw_paper(this_paper)
            n += 1

    temp += '<h2>' + str(n) + ' Publications From Unknown Years</h2>'
    temp += '<h3>Some data may be missing from these publications.</h3>'
    temp += html

    temp += build_common_foot()

    unknown_file.write(temp)


############################################################
# Build a list of all mesh keywords
############################################################
def build_mesh(papers):

    print("\n###HTML - mesh###")

    shutil.copyfile(config.template_dir + '/keyword_history.js', config.html_dir + '/mesh/keyword_history.js')

    mesh_papers_all = {}
    mesh_papers_major = {}
    html_file_all = open(config.html_dir + '/all_keywords/index.html', 'w')
    html_file_major = open(config.html_dir + '/major_keywords/index.html', 'w')

    # Build a dict of ALL mesh headings with a list of each hash (paper) that has this mesh term
    for this_paper in papers:
        try:
            # Look at all the mesh headings for this paper
            for this_mesh in this_paper['clean']['keywords']['mesh']:
                # If this mesh term is not already in the dict then add it
                if this_mesh['term'] not in mesh_papers_all:
                    mesh_papers_all[this_mesh['term']] = list()
                mesh_papers_all[this_mesh['term']].append(this_paper['IDs']['hash'])
        except:
            pass

    # Build a dict of ONLY MAJOR mesh headings with a list of each pmid in each
    data_from_count = 0
    for this_paper in papers:
        got_major_mesh = 0
        try:
            # Look at all the mesh headings for this paper
            for this_mesh in this_paper['clean']['keywords']['mesh']:
                # Only interested in major topics
                if this_mesh['major'] == 'Y':
                    got_major_mesh = 1
                    # If this mesh term is not in the dict then add it
                    if this_mesh['term'] not in mesh_papers_major:
                        mesh_papers_major[this_mesh['term']] = list()
                    mesh_papers_major[this_mesh['term']].append(this_paper['IDs']['hash'])

            data_from_count += got_major_mesh
        except:
            pass

    # # Read in mesh tree hierarchy
    # f = open(config.config_dir + "/mesh_tree_hierarchy.csv", 'rt')
    # mesh_tree = {}
    # mesh_tree_reverse = {}
    # try:
    #     reader = csv.reader(f)
    #     for row in reader:
    #         mesh_tree[row[2]] = row[0]
    #         mesh_tree_reverse[row[0]] = row[2]

    # finally:
    #     f.close()

    # f = open(config.config_dir + "/mesh_categories.csv", 'rt')
    # mesh_categories = {}
    # try:
    #     reader = csv.reader(f)
    #     for row in reader:
    #         mesh_categories[row[0]] = row[1]

    # finally:
    #     f.close()

    # # Make list of second level mesh headings
    # second_found = 0
    # top_found = 0
    # total = 0
    # mesh_second_level_headings = {}
    # mesh_top_level_headings = {}
    # for this_mesh in mesh_papers_major:
    #     try:
    #         tree_number = mesh_tree[this_mesh]
    #         # Get Second Level
    #         tree_number_split = tree_number.split(".")
    #         second_level = mesh_tree_reverse[tree_number_split[0]]

    #         try:
    #             mesh_second_level_headings[second_level] += len(mesh_papers_major[this_mesh])
    #         except:
    #             mesh_second_level_headings[second_level] = len(mesh_papers_major[this_mesh])
    #         second_found += len(mesh_papers_major[this_mesh])

    #         # Get Top Level
    #         top_level = tree_number[0]
    #         try:
    #             mesh_top_level_headings[mesh_categories[top_level]] += len(mesh_papers_major[this_mesh])
    #         except:
    #             mesh_top_level_headings[mesh_categories[top_level]] = len(mesh_papers_major[this_mesh])
    #         top_found += len(mesh_papers_major[this_mesh])
    #     except:
    #         pass

    #     total += len(mesh_papers_major[this_mesh])

    # print "MeSH Second Level Found: " + str(second_found) + "/" + str(total)
    # print "Unique MeSH Second Level: " + str(len(mesh_second_level_headings))
    # print "MeSH Top Level Found: " + str(top_found) + "/" + str(total)
    # print "Unique MeSH Top Level: " + str(len(mesh_top_level_headings))

    # print "\n" + str(mesh_top_level_headings)
    # print "\n" + str(mesh_second_level_headings)

    # print "Second Level MeSH"
    # for mesh in mesh_second_level_headings:
        # print mesh + "\t" + str(mesh_second_level_headings[mesh])

    # print "\n\n"

    # print "Top Level MeSH"
    # for mesh in mesh_top_level_headings:
        # print mesh + "\t" + str(mesh_top_level_headings[mesh])

    # # Print mesh_papers
    # # Make a JSON file for each mesh term, in it put all the PMIDs for this term
    # for this_mesh in mesh_papers_all:
    #     file_name = config.html_dir + '/mesh/all_' + this_mesh
    #     fo = open(file_name, 'w')
    #     fo.write(json.dumps(mesh_papers_all[this_mesh], indent=4))
    #     fo.close()

    # # Make a JSON file for each major mesh term, in it put all the PMIDs for this term
    # for this_mesh in mesh_papers_major:
    #     file_name = config.html_dir + '/mesh/major_' + this_mesh
    #     fo = open(file_name, 'w')
    #     fo.write(json.dumps(mesh_papers_major[this_mesh], indent=4))
    #     fo.close()

    ######################################
    # Make HTML index page for ALL MESH headings
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; All Keywords</p>', "../", "")

    temp += '<h1 id="pagetitle">All Keywords</h1>'
    temp += '<p>' + str(len(mesh_papers_all)) + ' Keywords</p>'

    html_file_all.write(temp)

    # Make a page with ALL the headings on it
    html_file_all.write('<ul>')
    for this_mesh in sorted(mesh_papers_all):
        temp = '<li><a href="../mesh/' + this_mesh.replace(" ", "%20") + '/index.html">' + this_mesh + '</a></li>'
        html_file_all.write(temp)
    html_file_all.write('</ul>')

    temp = build_common_foot()
    html_file_all.write(temp)
    ######################################


    ######################################
    # Make HTML index page for MAJOR MESH headings
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Major Keywords (MeSH)</p>', "../", "")

    temp += '<h1 id="pagetitle">Major Keywords (MeSH)</h1>'
    temp += '<p>' + str(len(mesh_papers_major)) + ' Keywords</p>'

    html_file_major.write(temp)

    # Make a page with the MAJOR headings on it
    html_file_major.write('<ul>')
    for this_mesh in sorted(mesh_papers_major):
        temp = '<li><a href="../mesh/' + this_mesh.replace(" ", "%20") + '/index.html">' + this_mesh + '</a></li>'
        html_file_major.write(temp)
    html_file_major.write('</ul>')

    temp = build_common_foot()
    html_file_major.write(temp)
    ######################################

    # word_cloud_list = "["
    # word_cloud_n = 0
    # word_cloud_max = 0

    # word_cloud_raw = ""

  

    ############################################
    # Make an HTML page for ALL MESH terms
    for this_mesh in mesh_papers_all:

        if not os.path.exists(config.html_dir + '/mesh/' + this_mesh):
            os.mkdir(config.html_dir + '/mesh/' + this_mesh)

            ############################################
            # Calculate keyword usage and citations over time.
            # Saved to a json file for google chart to use
 
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
                        this_year = this_paper['clean']['clean_date']['year']
                        # Make sure there is a dict item for this year
                        if this_year not in summary:
                            summary[this_year] = {'num_papers': 0, 'citations': 0}

                        # increment the number of citations by one
                        summary[this_year]['num_papers'] += 1

                        # add the citations for this paper to the year running total
                        try:
                            summary[this_year]['citations'] += int(this_paper['clean']['citations']['scopus']['count'])
                        except:
                            pass

                    except:
                        pass

                # Add in some zeros when there is no papers for this year
                years = list(summary.keys())
                first_year = min(years)
                last_year = max(years)
                for this_year in range(int(first_year), int(last_year)):
                    try:
                        summary[str(this_year)]['num_papers']
                    except:
                        summary[str(this_year)] = {'num_papers': 0, 'citations': 0}

            # Print data to file
            data_file = open(config.html_dir + '/mesh/' + this_mesh + '/stats.js', 'w')
            # print >>data_file, 'var papers =([[\'Year\', \'Number of papers\'],'
            data_file.write('var papers =([[\'Year\', \'Number of papers\'],')
            for this_year in sorted(summary, reverse=False):
                data_file.write('[\''+this_year+'\','+str(summary[this_year]['num_papers'])+'],')
            data_file.write(']);')

            data_file.write('var citations =([[\'Year\', \'Number of Citations\'],')
            for this_year in sorted(summary, reverse=False):
                data_file.write('[\''+this_year+'\','+str(summary[this_year]['citations'])+'],')
            data_file.write(']);')

        # Output the HTML for this mesh term
        file_name = config.html_dir + '/mesh/' + this_mesh + '/index.html'
        with codecs.open(file_name, 'wb', "utf-8") as fo:

            # Put html together for this page
            temp = '<!DOCTYPE html><html lang="en-GB">'

            # html head
            temp += '<head>'
            temp += '<title>' + site_second_title + '</title>'
            temp += '<link rel="stylesheet" href="../../css/style_main.css">'
            temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'

            temp += '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
            # temp += '<script type="text/javascript" src="../' + this_mesh.replace(" ", "%20") + '.js"></script>'
            temp += '<script type="text/javascript" src="stats.js"></script>'
            temp += '<script>var primary_colour = "#' + config.project_details['colour_hex_primary'] + '";</script>'
            temp += '<script>var secondary_colour = "#' + config.project_details['colour_hex_secondary'] + '";</script>'
            temp += '<script type="text/javascript" src="../keyword_history.js"></script>'

            temp += '</head>'

            temp += build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; Keyword &gt; ' + this_mesh + '</p>', "../../", "")

            temp += '<h1 id="pagetitle">Keyword - ' + this_mesh + '</h1>'
            temp += '<h2>Keyword History</h2>'

            # Placeholers for the charts
            temp += '<div id="papers_chart_div"></div>'
            temp += '<div id="citations_chart_div"></div>'

            # List publications
            temp += '<h2>Publications</h2>'
            temp += '<p><em>' + str(len(mesh_papers_all[this_mesh])) + ' publications with this keyword</em></p>'

            fo.write(temp)

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

                    # Call draw paper function
                    html = draw_paper(this_paper)

                    fo.write(html)

                except:
                    pass

            temp = build_common_foot()
            fo.write(temp)

        fo.close()


    # # MAJOR KEYWORD CLOUD
    # # Get first x MAJOR keywords in order
    # # Prepare variables
    # ordered_mesh_papers_all = {}
    # mesh_papers_all_temp = {}
    # for this_mesh in mesh_papers_major:
    #     mesh_papers_all_temp[this_mesh] = mesh_papers_all[this_mesh]

    # # Sort and get x words
    # for x in range(1, 150):
    #     max_val = -1
    #     max_mesh = None
    #     for this_mesh in mesh_papers_all_temp:
    #         if len(mesh_papers_all_temp[this_mesh]) > max_val:
    #             max_val = len(mesh_papers_all_temp[this_mesh])
    #             max_mesh = this_mesh

    #     if max_mesh is not None:
    #         ordered_mesh_papers_all[max_mesh] = mesh_papers_all_temp[max_mesh]
    #         mesh_papers_all_temp.pop(max_mesh, None)

    # # Make papers list for headings pages
    # for this_mesh in ordered_mesh_papers_all:

    #     # Word cloud
    #     if word_cloud_n > 0:
    #         word_cloud_list += ','
    #     word_cloud_n += 1

    #     number = len(mesh_papers_all[this_mesh])

    #     if number > word_cloud_max:
    #         word_cloud_max = number

    #     word_cloud_list += '{"text":"' + this_mesh + '", "size":' + str(math.sqrt(number)*2.5) + '}'
    #     # word_cloud_list += '{"text":"' + this_mesh + '", "size":' + str(number*5) + '}'

    #     for x in range(0, number):
    #         word_cloud_raw += " " + this_mesh

    # # word_cloud_list += "]"
    # # build_word_cloud(papers, word_cloud_list, data_from_count)


###########################################################
# Build a google map based on the lat longs provided before.
###########################################################
def build_google_map(papers):

    print("\n###HTML - Instititutions Map###")

    info = []
    number_of_points = 0
    for this_paper in papers:
        try:
            this_place = {'lat': this_paper['clean']['location']['latitude'], 'long': this_paper['clean']['location']['longitude'], 'name': this_paper['clean']['location']['clean_institute']}
            info.append(this_place)
            number_of_points += 1
        except:
            pass

    kml = "var locations =["
    for this_info in info:
        kml += '["' + this_info['name'] + '",' + str(this_info['lat']) + ',' + str(this_info['long']) + '],'
    kml += ']'

    kml_file = codecs.open(config.html_dir + '/map/map.kml', 'w', 'utf-8')
    kml_file.write(kml)

    html_file = open(config.html_dir + '/map/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="https://maps.googleapis.com/maps/api/js"></script>'
    temp += '<script type="text/javascript" src="map.kml"></script>'
    temp += '<script type="text/javascript" src="map.js"></script>'

    shutil.copyfile(config.template_dir + '/map.js', config.html_dir + '/map/map.js')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/css/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Institutions Map</p>', "../", "onload='initialize()'")

    temp += '<h1 id="pagetitle">Institutions Map</h1>'

    temp += "<div id='map_canvas'></div>"
    temp += "<p>Data from " + intWithCommas(number_of_points) + " publications</p>"

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
# Publications by country
###########################################################
def build_country_map(papers):

    print("\n###HTML - Country Map###")

    countries = {}
    number_of_points = 0
    for this_paper in papers:
        try:

            if this_paper['clean']['location']['country'] in countries:
                countries[this_paper['clean']['location']['country']] += 1
            else:
                countries[this_paper['clean']['location']['country']] = 1
            number_of_points += 1
        except:
            pass

    country_string = ""
    for country in list(countries.keys()):
        country_string += ',["' + country + '",' + str(countries[country]) + ']'

    # Google charts uses different names than wikidata for some countries.
    country_string = country_string.replace("People's Republic of China", "China")
    country_string = country_string.replace("United States of America", "United States")

    html_file = open(config.html_dir + '/country/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script> <script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    temp += '<script type="text/javascript">' + "google.charts.load('current', {'packages':['geochart']});google.charts.setOnLoadCallback(drawRegionsMap);function drawRegionsMap() {var data = google.visualization.arrayToDataTable([ ['Country', 'Publications']" + country_string + "]); var options = { colorAxis: {colors: ['#" + config.project_details['colour_hex_secondary'] + "', '#" + config.project_details['colour_hex_primary'] + "']} }; var chart = new google.visualization.GeoChart(document.getElementById('regions_div')); chart.draw(data, options); }</script>"

    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/country/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/country/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications by Country</p>', "../", "")

    temp += '<h1 id="pagetitle">Publications by Country</h1>'

    temp += '<div id="regions_div" style="width: 900px; height: 500px;"><img src="loading.gif" alt="Loading"></div>'
    temp += "<p>Data from " + intWithCommas(number_of_points) + " publications</p>"

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
# Publications by UK institute
###########################################################
def build_institute_map(papers):

    print("\n###HTML - Institute Map###")

    institutes = {}
    number_of_points = 0
    for this_paper in papers:
        try:
            if this_paper['clean']['location']['clean_institute'] != "" and this_paper['clean']['location']['latitude'] != "" and this_paper['clean']['location']['longitude'] != "":
                if this_paper['clean']['location']['clean_institute'] in institutes:
                    institutes[this_paper['clean']['location']['clean_institute']]['count'] += 1
                else:
                    institutes[this_paper['clean']['location']['clean_institute']] = {}
                    institutes[this_paper['clean']['location']['clean_institute']]['lat'] = this_paper['clean']['location']['latitude']
                    institutes[this_paper['clean']['location']['clean_institute']]['lon'] = this_paper['clean']['location']['longitude']
                    institutes[this_paper['clean']['location']['clean_institute']]['count'] = 1
                number_of_points += 1
        except:
            pass

    institute_string = ""
    for this_institute in list(institutes.keys()):
        # note that we encode as ascii with 'ignore' set so all non-ascii
        # chars are removed. this is then decoded back into utf-8
        institute_string += ',[' + institutes[this_institute]['lat'] + ',' + institutes[this_institute]['lon'] + ',"' + str(this_institute.encode('ascii', 'ignore').decode('ascii')) + '",' + str(institutes[this_institute]['count']) + ']'

    html_file = open(config.html_dir + '/institute/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>'
    temp += "<script>google.charts.load('current', {'packages':['geochart']});google.charts.setOnLoadCallback(drawMarkersMap);function drawMarkersMap() {var data = google.visualization.arrayToDataTable([['lat', 'lon', 'Institute','Publication count']" + institute_string + " ]); var options = {magnifyingGlass: {zoomFactor: '15.0'}, region: 'GB', displayMode: 'markers', colorAxis: {colors: ['#" + config.project_details['colour_hex_secondary'] + "', '#" + config.project_details['colour_hex_primary'] + "']}}; var chart = new google.visualization.GeoChart(document.getElementById('regions_div'));chart.draw(data, options); };</script>"

    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/institute/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/institute/map.css')

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications by UK City</p>', "../", "")

    temp += '<h1 id="pagetitle">Publications by UK Institute</h1>'

    temp += '<div id="regions_div" style="width: 900px; height: 500px;"><img src="loading.gif" alt="Loading"></div>'
    temp += "<p>Data from " + intWithCommas(number_of_points) + " publications</p>"
    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
# Util to stick in commas between 1000s in big integers
###########################################################
def intWithCommas(x):
    if not isinstance(x, int):
        raise TypeError("Parameter must be an integer.")
    if x < 0:
        return '-' + intWithCommas(-x)
    result = ''
    while x >= 1000:
        x, r = divmod(x, 1000)
        result = ",%03d%s" % (r, result)
    return "%d%s" % (x, result)


###########################################################
# Build metrics page
###########################################################
def build_metrics(papers, age_weighted_citation, age_weighted_citation_data, study_start_year, study_current_year):

    print("\n###HTML - Metrics###")

    html_file = open(config.html_dir + '/metrics/index.html', 'w')

    # Metric calculations
    total_publications = len(papers)
    total_citations = 0
    total_citations_data_from_count = 0
    paper_citations = []
    c20_index = 0

    study_duration = study_current_year - study_start_year

    c_index_bound = 100

    for this_paper in papers:
        # add the citations for paper to total
        try:
            cit = int(this_paper['clean']['citations']['scopus']['count'])
            total_citations += cit
            total_citations_data_from_count += 1
            paper_citations.append(cit)

            # increment c20-index if more that 20 citations
            if cit >= c_index_bound:
                c20_index += 1
        except:
            pass

    # We might not have any citattions - e.g. if quotas hit.
    if total_citations > 0 and total_citations_data_from_count > 0:
        average_citations = float(total_citations)/float(total_citations_data_from_count)
    else:
        average_citations = 0

    i20_index_per_year = float(c20_index)/float(study_duration)

    # cal h-index
    paper_citations.sort(reverse=True)
    h_index = 0
    # cits_so_far = 0

    for x in range(0, len(paper_citations)):
        if x > paper_citations[x]:
            break
        h_index = x
        # cits_so_far += paper_citations[x]

    # cal g-index
    g_index = 0
    cits_so_far = 0

    for x in range(0, len(paper_citations)):
        cits_so_far += paper_citations[x]
        if cits_so_far < x * x:
            break
        g_index = x

    # NUMBER OF PAPERS PER CITATION COUNT
    citation_number_limit = 200
    citation_bin_size = 20

    num_papers_citations = []
    max_citations = 0
    list_of_citation_counts = []

    # Get the max number of citations
    for this_paper in papers:
        try:
            n_cits = int(this_paper['clean']['citations']['scopus']['count'])
            if n_cits > max_citations:
                max_citations = n_cits
        except:
            pass

    # Create a array of zeros of length max_citations + 1
    for x in range(0, max_citations + 1):
        num_papers_citations.insert(x, 0)

    # Count citations into array
    for this_paper in papers:
        try:
            n_cits = int(this_paper['clean']['citations']['scopus']['count'])
            list_of_citation_counts.append(n_cits)
            num_papers_citations[n_cits] += 1
        except:
            pass

    # Get the median number of citations
    if len(list_of_citation_counts) > 0:
        list_of_citation_counts.sort()
        median_citations = list_of_citation_counts[round(len(list_of_citation_counts)/2)]
    else:
        median_citations = 0

    # = Low Citations Range =
    # Create data string for plot
    n_papers_with_x_citations = "var papers_per_citation_count = ([['Number of Citations (Scopus)','Number of Papers',{ role: 'style' }]"
    for this_n_citations in range(1, citation_number_limit):

        # colour used to mark average value
        colour = ""
        if this_n_citations == round(average_citations, 0):
            colour = "#" + config.project_details['colour_hex_secondary']
        if this_n_citations == median_citations:
            colour = "green"

        try:
            n_papers_with_x_citations += ",[" + str(this_n_citations) + "," + str(num_papers_citations[this_n_citations]) + ",'" + colour + "']"
        except:
            n_papers_with_x_citations += ",[" + str(this_n_citations) + ",0,'" + colour + "']"

    n_papers_with_x_citations += "]);"

    # High Citations Range, but only if they are above the citation_number_limit
    if(max_citations >= citation_number_limit):
        n_papers_with_x_citations += "var papers_per_high_citation_count = ([['Number of Citations (Scopus)','Number of Papers',{ role: 'style' }]"
        for this_bin in range(0, round((max_citations - citation_number_limit)/citation_bin_size) + 1):

            # Calculate the start and end of the bin
            bin_start = this_bin * citation_bin_size + citation_number_limit + 1
            bin_end = bin_start + citation_bin_size

            num_papers_in_bin = 0

            # Count all citation counts in bin
            for n in range(bin_start, bin_end):
                try:
                    num_papers_in_bin += num_papers_citations[n]
                except:
                    pass

            n_papers_with_x_citations += ",['" + str(bin_start) + "-" + str(bin_end) + "'," + str(num_papers_in_bin) + ",'" + colour + "']"

        n_papers_with_x_citations += "]);"

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
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

    temp += '(All citation calculations based on <a href="https://www.scopus.com">Scopus</a> data.)'

    # Output Metrics
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
    temp += "<div class='metric_name'>Mean Citations Per Publication</div>"
    temp += "<div class='metric_value'>" + str("{0:.2f}".format(round(average_citations, 2))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The total number of citations divided by the total number of publications.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Age-weighted Mean Citations Per Publication</div>"
    temp += "<div class='metric_value'>" + str("{0:.2f}".format(round(age_weighted_citation, 3))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + intWithCommas(age_weighted_citation_data) + " Publications</div>"
    temp += "<div class='metric_description'>Age-weighted Mean Citations Per Publication.</div>"
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
    temp += "<p style='text-align:center;'>Data from " + intWithCommas(age_weighted_citation_data) + " publications</p>"
    temp += '<div id="papers_per_year_div"></div>'
    temp += "<p style='text-align:center;'>Data from " + intWithCommas(age_weighted_citation_data) + " publications</p>"
    temp += '<div id="papers_per_citation_count_div"></div>'
    temp += "<div style='margin-left:auto;margin-right:auto;'><div class='average_citations' style='height:15px; width:33px; float:left; background:#" + config.project_details['colour_hex_secondary'] + "'></div><div style='height: 15px;line-height: 15px;padding-left: 40px;'> Mean number of citations</div></div>"
    temp += "<div style='margin-left:auto;margin-right:auto;margin-top:5px;'><div class='average_citations' style='height:15px; width:33px; float:left; background:green'></div><div style='height: 15px;line-height: 15px;padding-left: 40px;'> Median number of citations</div></div>"
    temp += "<p style='text-align:center;'>Data from " + intWithCommas(total_citations_data_from_count) + " publications</p>"

    if(max_citations >= citation_number_limit):
        temp += '<div id="papers_per_high_citation_count_div"></div>'
        temp += "<p style='text-align:center;'>Data from " + intWithCommas(total_citations_data_from_count) + " publications</p>"

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


# ###########################################################
# # Build keyword word cloud
# ###########################################################
# def build_word_cloud(papers, list, data_from_count):

#     print("\n###HTML - Keyword Cloud###")

#     html_file = open(config.html_dir + '/wordcloud/index.html', 'w')

#     # Put html together for this page
#     temp = '<!DOCTYPE html><html lang="en-GB">'

#     # html head
#     temp += '<head>'
#     temp += '<title>' + site_second_title + '</title>'
#     temp += '<link rel="stylesheet" href="../css/style_main.css">'
#     temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
#     temp += '<link rel="stylesheet" href="../css/map.css">'

#     shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/wordcloud/d3wordcloud.js')
#     shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/wordcloud/d3.layout.cloud.js')

#     temp += '<script>var word_list = ' + list + '</script>'
#     temp += '<script src="https://d3js.org/d3.v3.min.js"></script>'
#     temp += '<script src="d3.layout.cloud.js"></script>'

#     temp += '</head>'

#     temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Major Keyword Cloud</p>', "../", "")

#     temp += '<h1 id="pagetitle">Major Keyword Cloud</h1>'

#     temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
#     temp += '</cloud>'
#     temp += "<p>Data from " + intWithCommas(data_from_count) + " publications</p>"

#     temp += '<script src="d3wordcloud.js"></script>'

#     html_file.write(temp)

#     temp = build_common_foot()
#     html_file.write(temp)


###########################################################
# Build abstract word cloud
###########################################################
def build_abstract_word_cloud(papers, data_from_count):

    print("\n###HTML - Abstract Word Cloud###")

    biggest_word_size = 100
    abstract_words = {}
    d3_word_list = "["
    n = 0

    with open(config.data_dir + "/abstract_raw.csv", 'rt') as f:

        reader = csv.reader(f)
        for row in reader:
            try:
                if row[0] != "":
                    abstract_words[row[0]] = int(row[1])
            except:
                pass

    # Sort the words
    abstract_words = sorted(abstract_words.items(), key=lambda x: x[1], reverse=True)

    # Only pick the top 200
    abstract_words = abstract_words[:200]

    # Grab the most frequent word
    max_count = abstract_words[0][1]

    # Scale the words so the biggest one is equal to biggest_word_size
    abstract_words = {k: round(v * biggest_word_size / max_count) for k, v in abstract_words}

    for this_word in abstract_words:
        if n > 0:
            d3_word_list += ","

        # word_list += '["' + row[0].replace("'","\'").replace('"','\"') + '",' + str(row[1]) +  ']'
        # word_list += '{"text":"' + row[0].replace("'", "\'").replace('"', '\"') + '","size":' + str(math.sqrt(int(row[1]))*1.5) + '}'
        # word_list += '{"text":"' + row[0].replace("'", "\'").replace('"', '\"') + '","size":' + str(row[1]) + '}'

        d3_word_list += '{"text":"' + this_word.replace("'", "\'").replace('"', '\"') + '","size":' + str(abstract_words[this_word]) + '}'
        n += 1

    d3_word_list += "];"
    list_file = open(config.html_dir + '/abstractwordcloud/list.js', 'w')
    list_file.write(" var word_list = " + d3_word_list)

    html_file = open(config.html_dir + '/abstractwordcloud/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'
    # temp += '<style>.wordcloud{ width:100%; height:500px;}</style>'

    shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/abstractwordcloud/d3wordcloud.js')
    shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/abstractwordcloud/d3.layout.cloud.js')

    temp += '<script src="list.js"></script>'
    temp += '<script src="https://d3js.org/d3.v3.min.js"></script>'
    temp += '<script src="d3.layout.cloud.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Abstract Word Cloud</p>', "../", "")

    temp += '<h1 id="pagetitle">Abstract Word Cloud</h1>'

    temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
    temp += '</cloud>'

    temp += '<script src="d3wordcloud.js"></script>'
    temp += "<p>Data from " + intWithCommas(data_from_count) + " publications</p>"

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)



###########################################################
# Build keyword word cloud
###########################################################
def build_keyword_word_cloud(papers, data_from_count):

    print("\n###HTML - Abstract Word Cloud###")

    biggest_word_size = 100
    words = {}
    d3_word_list = "["
    n = 0

    with open(config.data_dir + "/keywords_raw.csv", 'rt') as f:
        reader = csv.reader(f)
        for row in reader:
            try:
                if row[0] != "":
                    words[row[0]] = int(row[1])
            except:
                pass

    # Sort the words
    words = sorted(words.items(), key=lambda x: x[1], reverse=True)

    # Only pick the top 200
    words = words[:200]

    # Grab the most frequent word
    max_count = words[0][1]

    # Scale the words so the biggest one is equal to biggest_word_size
    words = {k: round(v * biggest_word_size / max_count) for k, v in words}

    for this_word in words:
        if n > 0:
            d3_word_list += ","

        d3_word_list += '{"text":"' + this_word.replace("'", "\'").replace('"', '\"') + '","size":' + str(words[this_word]) + '}'
        n += 1

    d3_word_list += "];"
    list_file = open(config.html_dir + '/keyword_wordcloud/list.js', 'w')
    list_file.write(" var word_list = " + d3_word_list)

    html_file = open(config.html_dir + '/keyword_wordcloud/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'
    
    shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/keyword_wordcloud/d3wordcloud.js')
    shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/keyword_wordcloud/d3.layout.cloud.js')

    temp += '<script src="list.js"></script>'
    temp += '<script src="https://d3js.org/d3.v3.min.js"></script>'
    temp += '<script src="d3.layout.cloud.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Keyword Word Cloud</p>', "../", "")

    temp += '<h1 id="pagetitle">Keyword Word Cloud</h1>'

    temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
    temp += '</cloud>'

    temp += '<script src="d3wordcloud.js"></script>'
    temp += "<p>Data from " + intWithCommas(data_from_count) + " publications</p>"

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
#
###########################################################
def get_author_string_from_hash(hash_string, network):

    for author in network['authors']:
        if author == hash_string:
            return network['authors'][author]['clean']


###########################################################
# Build Author Network
###########################################################
def build_author_network(papers, network):

    print("\n###HTML - Author Network###")

    # Create json file
    net_file = open(config.html_dir + '/authornetwork/network.json', 'w')

    net_json = '{'
    net_json += '"nodes":['
    n = 0
    net_file.write(net_json)

    for author in network['authors']:
        net_json = ""
        if n > 0:
            net_json += ","
        net_json += '{"id": "' + network['authors'][author]['clean'] + '", "group":1}'

        net_file.write(net_json)
        n += 1

    net_json = '],"links": ['
    net_file.write(net_json)

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

            net_file.write(net_json)
            n += 1
        except:
            pass

    net_json = ']'
    net_json += '}'
    net_file.write(net_json)

    html_file = open(config.html_dir + '/authornetwork/index.html', 'w')

    shutil.copyfile(config.template_dir + '/network.js', config.html_dir + '/authornetwork/network.js')

    if config.page_show_author_network:
        try:
            shutil.copyfile(config.config_dir + '/' + config.project_details['short_name'] + '_author_network.png', config.html_dir + '/authornetwork/author_network.png')
        except:
            logging.warn("Not Author Network Image")

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<style>.links line {  stroke: #999;  stroke-opacity: 0.6;} .nodes circle {  stroke: #fff;  stroke-width: 1.5px;} </style>'

    temp += '<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.1.0/jquery.min.js"></script>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Author Network</p>', "../", "")

    temp += '<h1 id="pagetitle">Author Network</h1>'

    # Print nodes to csv
    nodes_csv = open(config.html_dir + '/authornetwork/nodes.csv', 'w')

    nodes_csv.write('id,Label')
    n = 0

    for author in network['authors']:
        nodes_csv.write(author + "," + network['authors'][author]['clean'])
        n += 1

    # Print connections to csv
    connections_csv = open(config.html_dir + '/authornetwork/connections.csv', 'w')

    connections_csv.write('Source,Target')

    n = 0
    for con in network['connections']:
        try:

            author_0 = network['connections'][con]['authors'][0]['author_hash']
            author_1 = network['connections'][con]['authors'][1]['author_hash']

            n_con = network['connections'][con]['num_connections']/2

            connections_csv.write('"' + author_0 + '","' + author_1 + '"')

        except:
            pass
        n += 1

    temp += '<a id="network" href="author_network.png"><img src="author_network.png" alt="Author Network"></a>'
    temp += '<p style="display:none;" id="no_network">No Author Network Image.</p>'
    temp += "<script>var xmlhttp = new XMLHttpRequest();xmlhttp.onreadystatechange = function() {if (xmlhttp.readyState == 4 && xmlhttp.status == 404) {document.getElementById('network').style.display = 'none';document.getElementById('no_network').style.display = 'block';}};xmlhttp.open('GET', 'author_network.png', true);xmlhttp.send();</script>"

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
# Help Page
###########################################################
def build_help():

    print("\n###HTML - Help/Information page###")

    html_file = open(config.html_dir + '/help/index.html', 'w')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    # temp += '<link rel="stylesheet" href="../css/map.css">'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Information</p>', "../", "")

    temp += '<h1 id="pagetitle">Information</h1>'

    temp += '<h2>Where does the data come from?</h2>'
    temp += '<h3>Publication Data</h3>'
    temp += '<p>Data and metadata for the publications are downloaded using the PubMed and DOI APIs.</p>'

    temp += '<h3>What does PubMed, DOI and API mean?</h3>'
    temp += '<p>PubMed is an online search engine which contains journal articles on lifesciences and biomedical topics. DOI is short for Digital Object Identifier and is the de facto way of referring to journal articles, they look like  <a href="https://doi.org/10.12688/f1000research.25484.2">https://doi.org/10.12688/f1000research.25484.2</a> and <a href="http://doi.org">http://doi.org</a> are in charge of them. An API is an Application Programming Interface, which is a fancy way of saying that computer programs can send and receive data to a service. So here we send and receive data to/from PubMed and doi.org.</p>'

    temp += '<h3>Citations Data</h3>'
    temp += '<p>Citation data is retrieved from the <a href="https://www.elsevier.com/solutions/scopus">Scopus</a> API provided by Elsevier and from the <a href="https://europepmc.org/">Europe PMC</a> API.</p>'

    temp += '<h3>Geodata</h3>'
    temp += '<p>The coordinate location and town of institutions is retrieved from the free project <a href="https://www.wikidata.org">Wikidata</a>.</p>'

    temp += '<h2>What are the circles with numbers in them?</h2>'
    temp += '<p>Journal article citations are one way to track how an article is being used. Another way is by considering a wider set of metrics, such as if it is mentioned in news articles, or on social media. This is what <a href"https://www.altmetric.com">https://www.altmetric.com</a> does, and by hovering your mouse over one of the circles (or clicking on it) you can get an overview of where this article is being talked about.</p>'

    temp += "<h2>Why don't some statistics use data from all publications?</h2>"

    temp += '<p>Data can be missing because some publications are very old and the data about it just doesn\'t exist, or sometimes a recent paper may not have the data about it available <em>yet</em>. There are many different ways to track publications'
    temp += ' and prior to the introduction of DOIs in 2000 there was no standard method. This means some metadata on old publications could have been lost or not recorded.</p>'

    temp += '<p>The data used for the statistics are gathered from databases which only collect data from particular journals.'
    temp += ' This problem is most obvious for citation data from Scopus and Europe PMC. Although they have a particular publication on record,'
    temp += ' there may be citations for this publication from a journal that they do not index and therefore these will not be in the citation count. This is why citation counts from different sources are not always the same.</p>'

    temp += '<h2>What\' the difference between "All keywords" and "Major keywords (MeSH)"?</h2>'
    temp += '<p>These keywords are assigned to journal articles when they are published. MeSH (Medical Subject Headings) is a commonly used hierarchical controlled vocabulary. This means that there is a list of well defined terms which are allowed to be used (controlled vocabulary), and they are related to one another (hierarchical) e.g. "finger" belongs to "hand". So the Major Keywords are the broader areas, and all keywords will be more fine grained. </p>'



    temp += '<h2>I want to know more!</h2>'
    temp += '<ul>'
    temp += '<li>The source code is at <a href="https://github.com/OllyButters/puma">https://github.com/OllyButters/puma</a></li>'
    temp += '<li>Some documentation is at <a href="https://github.com/OllyButters/puma/wiki">https://github.com/OllyButters/puma/wiki</a></li>'
    temp += '<li>There is even a paper written about it: <a href="https://f1000research.com/articles/9-1095/v2">https://f1000research.com/articles/9-1095/v2</a></li>'
    temp += '<li>You can talk to us at: <a href="https://twitter.com/DrOllyButters">@DrOllyButters</a>, <a href="https://twitter.com/DrBeccaWilson">@DrBeccaWilson</a> and <a href="https://twitter.com/_hugh_garner_">@_hugh_garner_</a></li></ul>'
    
    temp += '<h2>Funding</h2>'
    temp += 'This project has been funded by:'
    temp += '<ul>'
    temp += '<li>CLOSER, whose mission is to maximise the use, value and impact of longitudinal studies. CLOSER is funded by the Economic and Social Research Council (ESRC) and Medical Research Council (MRC) (grant reference: ES/K000357/1).</li>'
    temp += '<li>Becca Wilson is a UKRI Innovation Fellow with HDR UK [MR/S003959/1].</li>'
    temp += '<li>The Nuffield Foundation research placement program.</li>'
    temp += '<li>The Wellcome Trust and Medical Research Council (grant number 108439/Z/15/Z).</li>'
    temp += '<li>The European Union’s Horizon 2020 research and innovation programme under grant agreement No 824989 and the Canadian Institutes of Health Research (CIHR).</li>'
    temp += '<li>The National Institute for Health Research Applied Research Collaboration.</li>'
    temp += '</ul>'

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
# Search Page
###########################################################
def build_search(papers):

    print("\n###HTML - Search page###")

    html_file = open(config.html_dir + '/search/index.html', 'w')
    shutil.copyfile(config.template_dir + '/search.js', config.html_dir + '/search/search.js')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + site_second_title + '</title>'
    temp += '<link rel="stylesheet" href="../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../css/colour_scheme.css">'
    temp += '<link rel="stylesheet" href="../css/map.css">'
    temp += '<script>var primary_colour = "' + config.project_details['colour_hex_primary'] + '";</script>'
    temp += '<script>var secondary_colour = "' + config.project_details['colour_hex_secondary'] + '";</script>'
    temp += '<script>var name = "' + config.project_details['name'] + '";</script>'
    temp += '<script src="search.js"></script>'

    temp += '<style> button { height:2em; font-size:1em; margin-left:5px; } input#search { width:50%; }</style>'

    temp += '</head>'

    temp += build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Search</p>', "../", "")

    temp += '<h1 id="pagetitle">Search</h1>'

    temp += '<p>Search for the fields titles, abstracts, keywords and authors.<br/>To narrow down the results search for multiple fields at once.</p>'
    temp += '<p><input type="text" id="search"><button onclick="search();">Search</button></p>'

    temp += '<div style="display:none;" id="searching">Searching...</div>'
    temp += '<h2 id="num_search_results"></h2>'
    temp += '<div id="search_results"></div>'

    searchable_data = []
    for this_paper in papers:
        this_subset = {}
        this_subset['IDs'] = this_paper['IDs']
        this_subset['clean'] = {}
        try:
            this_subset['clean']['title'] = this_paper['clean']['title']
        except:
            pass

        try:
            this_subset['clean']['abstract'] = this_paper['clean']['abstract']
        except:
            pass

        try:
            this_subset['clean']['keywords'] = this_paper['clean']['keywords']
        except:
            pass

        try:
            this_subset['clean']['full_author_list'] = []
            for this_author in this_paper['clean']['full_author_list']:
                this_subset['clean']['full_author_list'].append({'clean': this_author['clean']})
        except:
            pass

        try:
            this_subset['clean']['journal'] = this_paper['clean']['journal']
        except:
            pass

        searchable_data.append(this_subset)

    temp += '<script id="search_data">var papers = ' + str(json.dumps(searchable_data)).replace("<", "&lt;").replace(">", "&gt;") + ';</script>'

    html_file.write(temp)

    temp = build_common_foot()
    html_file.write(temp)


###########################################################
# CSS colour scheme
###########################################################
def build_css_colour_scheme():

    # This function generates the CSS used for the colour scheme for the whole site.
    # This is includes the images for the top bar and naviagation bar.
    # The colours and image data is taken from the config file.

    html_file = open(config.html_dir + '/css/colour_scheme.css', 'w')

    temp = ".header-container {background: #" + config.project_details['colour_hex_primary'] + ";}"

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

    temp += ".header-container:before,"
    temp += ".header-container:after {"
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
    temp += "div#content ul a {"
    temp += "color: #" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    temp += "/* Metrics Page */"

    temp += "div#content div.metric {"
    temp += "background:#efede9;"
    temp += "color:#" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    # Hide the description text by default
    temp += "div#content div.metric div.metric_description {display: none;}"

    temp += "div#content div.metric:hover {"
    temp += "background:#" + config.project_details['colour_hex_primary'] + ";"
    temp += "}"

    temp += "div#content div.metric:hover div.metric_name {"
    temp += "color:#f2f2f2;"
    temp += "}"

    temp += "div#content div.metric:hover div.metric_value {"
    temp += "color:#f2f2f2;"
    temp += "}"

    temp += "div#content div.metric:hover div.metric_description {"
    temp += "color:#efede9; display: block;"
    temp += "}"

    temp += "a:link {color:#" + config.project_details['colour_hex_primary'] + "}"
    temp += "a:visited {color:#" + config.project_details['colour_hex_primary'] + "}"

    html_file.write(temp)
