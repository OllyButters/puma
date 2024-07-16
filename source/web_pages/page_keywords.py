import shutil
import os
import codecs

from config import config
from . import common_html as ch



############################################################
# Build a list of all mesh keywords and make a page for each
############################################################
def build_mesh(papers):

    print("\n###HTML - mesh###")

    mesh_papers_all = {}
    mesh_papers_major = {}
    other_keywords = {}
    html_file_all = open(config.html_dir + '/keywords/index.html', 'w', encoding='utf-8')
    html_file_major = open(config.html_dir + '/mesh/index.html', 'w', encoding='utf-8')

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

    # Build a dict of OTHER keywords with a list of each hash (paper) that has this keyword
    for this_paper in papers:
        try:
            # Look at all the OTHER keywords for this paper
            for this_keyword in this_paper['clean']['keywords']['other']:
                # If this keyword is not already in the dict then add it
                if this_keyword not in other_keywords:
                    other_keywords[this_keyword] = list()
                other_keywords[this_keyword].append(this_paper['IDs']['hash'])
        except:
            pass

    # Merge the MESH and other keywords together. Note that there are a load of
    # formatting issues to deal with, e.g. other keywords seem to be lower case.
    all_keywords = {}
    for this_mesh in mesh_papers_all:

        formatted_mesh = this_mesh[0].capitalize() + this_mesh[1:]

        all_keywords[formatted_mesh] = mesh_papers_all[this_mesh]

    for this_keyword in other_keywords:

        formatted_keyword = this_keyword[0].capitalize() + this_keyword[1:]

        if formatted_keyword in all_keywords.keys():
            all_keywords[formatted_keyword].extend(other_keywords[this_keyword])

        else:
            all_keywords[formatted_keyword] = other_keywords[this_keyword]

    mesh_papers_all = all_keywords

    ######################################
    # Make HTML index page for ALL keywords

    html = ch.build_common_head("../", "")
    html += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; All Keywords</p>', "../")

    html += '<h1 id="pagetitle">All Keywords</h1>'
    html += '<p>' + str(len(all_keywords)) + ' Keywords</p>'

    html_file_all.write(html)

    # Make a page with ALL the headings on it
    html_file_all.write('<ul>')
    for this_keyword in sorted(all_keywords):
        html = '<li><a href="../keywords/' + this_keyword.replace(" ", "%20").replace("/", "-") + '/index.html">' + this_keyword + '</a></li>'
        html_file_all.write(html)
    html_file_all.write('</ul>')

    html = ch.build_common_foot("../")
    html_file_all.write(html)
    ######################################

    ######################################
    # Make HTML index page for MAJOR MESH headings

    html = ch.build_common_head("../", "")
    html += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Major Keywords (MeSH)</p>', "../")

    html += '<h1 id="pagetitle">Major Keywords (MeSH)</h1>'
    html += '<p>' + str(len(mesh_papers_major)) + ' Keywords</p>'

    html_file_major.write(html)

    # Make a page with the MAJOR headings on it
    html_file_major.write('<ul>')
    for this_mesh in sorted(mesh_papers_major):
        html = '<li><a href="../mesh/' + this_mesh.replace(" ", "%20").replace("/", "-") + '/index.html">' + this_mesh + '</a></li>'
        html_file_major.write(html)
    html_file_major.write('</ul>')

    html = ch.build_common_foot("../")
    html_file_major.write(html)
    ######################################

    _make_keywords_pages(papers, mesh_papers_major, "mesh")
    _make_keywords_pages(papers, all_keywords, "keywords")


# Helper function to build the list of papers that have this keyword. Note that
# a mesh heading is consider a keyword here.
def _make_keywords_pages(papers, keywords, url_part):

    shutil.copyfile(config.template_dir + '/keyword_history.js', config.html_dir + '/' + url_part + '/keyword_history.js')

    ############################################
    # Make an HTML page for ALL MESH terms
    for this_keyword in keywords:


        # Some keywords have a / in them, which messes with URL and file paths. Just swap it out for a - here.        
        if this_keyword.find("/") > 0:
            this_keyword_safe = this_keyword.replace("/", "-")
        else:
            this_keyword_safe = this_keyword


        if not os.path.exists(config.html_dir + '/' + url_part + '/' + this_keyword_safe):
            os.mkdir(config.html_dir + '/' + url_part + '/' + this_keyword_safe)

            ############################################
            # Calculate keyword usage and citations over time.
            # Saved to a json file for google chart to use

            summary = {}
            # Calculate the number of papers for each year
            for this_paper in keywords[this_keyword]:

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
                for this_year in range(int(first_year), int(last_year) + 1):
                    try:
                        summary[str(this_year)]['num_papers']
                    except:
                        summary[str(this_year)] = {'num_papers': 0, 'citations': 0}

            # Print data to file
            data_file = open(config.html_dir + '/' + url_part + '/' + this_keyword_safe + '/stats.js', 'w', encoding='utf-8')
            
            data_file.write('var papers =([[\'Year\', \'Number of papers\'],')
            for this_year in sorted(summary, reverse=False):
                data_file.write('[\''+this_year+'\','+str(summary[this_year]['num_papers'])+'],')
            data_file.write(']);')

            data_file.write('var citations =([[\'Year\', \'Number of Citations\'],')
            for this_year in sorted(summary, reverse=False):
                data_file.write('[\''+this_year+'\','+str(summary[this_year]['citations'])+'],')
            data_file.write(']);')

        # Output the HTML for this mesh term
        file_name = config.html_dir + '/' + url_part + '/' + this_keyword_safe + '/index.html'
        with codecs.open(file_name, 'wb', "utf-8") as fo:

            # Put html together for this page

            # extra JS needed for the plots of keywords
            extra_head = '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
            extra_head += '<script type="text/javascript" src="stats.js"></script>'
            extra_head += '<script>var title="keyword";</script>'
            extra_head += '<script>var primary_colour = "#' + config.project_details['colour_hex_primary'] + '";</script>'
            extra_head += '<script>var secondary_colour = "#' + config.project_details['colour_hex_secondary'] + '";</script>'
            extra_head += '<script type="text/javascript" src="../keyword_history.js"></script>'
            extra_head += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"


            temp = ch.build_common_head("../../", extra_head)
            temp += ch.build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; Keyword &gt; ' + this_keyword + '</p>', "../../")

            temp += '<h1 id="pagetitle">Keyword - ' + this_keyword + '</h1>'
            temp += '<h2>Keyword History</h2>'

            # Placeholers for the charts
            temp += '<div id="papers_chart_div"></div>'
            temp += '<div id="citations_chart_div"></div>'

            # List publications
            temp += '<h2>Publications</h2>'
            temp += '<p><em>' + str(len(keywords[this_keyword])) + ' publications with this keyword</em></p>'

            fo.write(temp)

            # Build the text needed for each paper
            for this_paper in keywords[this_keyword]:

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
                    html = ch.draw_paper(this_paper, "../../")

                    fo.write(html)

                except:
                    pass

            temp = ch.build_common_foot("../../")
            fo.write(temp)

        fo.close()

