import shutil
import os
import codecs

from config import config
from . import utils
from . import common_html as ch

############################################################
# Build a list of all zotero tags, and make a page for each
############################################################
def build_zotero_tags(papers):

    print("\n###HTML - Zotero tags###")

    shutil.copyfile(config.template_dir + '/keyword_history.js', config.html_dir + '/tags/keyword_history.js')

    zotero_tags = {}
    zotero_tags_counts = {}
    html_file = open(config.html_dir + '/tags/index.html', 'w', encoding='utf-8')

    # Build a dict of ALL zotero tags with a list of each hash (paper) that has this tag.
    for this_paper in papers:
        try:
            # Look at all the zotero tags for this paper
            for this_tag in this_paper['clean']['zotero_tags']:
                # If this tag is not already in the dict then add it
                if this_tag['tag'] not in zotero_tags:
                    zotero_tags[this_tag['tag']] = list()
                zotero_tags[this_tag['tag']].append(this_paper['IDs']['hash'])
        except:
            pass



    ############################################
    # Make an HTML page for all tags
    for this_tag in zotero_tags:

        zotero_tags_counts[this_tag] = dict()

        # Some tags have a / in them, which messes with URL and file paths. Just swap it out for a - here.
        if this_tag.find("/") > 0:
            this_tag_safe = this_tag.replace("/", "-")
        else:
            this_tag_safe = this_tag

        if not os.path.exists(config.html_dir + '/tags/' + this_tag_safe):
            os.mkdir(config.html_dir + '/tags/' + this_tag_safe)

            ############################################
            # Calculate keyword usage and citations over time.
            # Saved to a json file for google chart to use

            summary = {}
            # Calculate the number of papers for each year
            for this_paper in zotero_tags[this_tag]:

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

                        # increment the number of papers by one
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

            # Save this data so we can use it for the summary tag page
            total_papers = 0
            total_citations = 0
            for this_year in range(int(first_year), int(last_year) + 1):
                total_papers = total_papers + int(summary[str(this_year)]['num_papers'])
                total_citations = total_citations + int(summary[str(this_year)]['citations'])
            zotero_tags_counts[this_tag]['total_papers'] = total_papers
            zotero_tags_counts[this_tag]['total_citations'] = total_citations

            # Print data to file
            data_file = open(config.html_dir + '/tags/' + this_tag_safe + '/stats.js', 'w', encoding='utf-8')

            data_file.write('var papers =([[\'Year\', \'Number of papers\'],')
            for this_year in sorted(summary, reverse=False):
                data_file.write('[\''+this_year+'\','+str(summary[this_year]['num_papers'])+'],')
            data_file.write(']);')

            data_file.write('var citations =([[\'Year\', \'Number of Citations\'],')
            for this_year in sorted(summary, reverse=False):
                data_file.write('[\''+this_year+'\','+str(summary[this_year]['citations'])+'],')
            data_file.write(']);')

        # Output the HTML for this tag
        file_name = config.html_dir + '/tags/' + this_tag_safe + '/index.html'
        with codecs.open(file_name, 'wb', "utf-8") as fo:

            # Put html together for this page

            # extra JS needed for the plots of tags
            extra_head = '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
            extra_head += '<script type="text/javascript" src="stats.js"></script>'
            extra_head += '<script>var title="tag";</script>'
            extra_head += '<script>var primary_colour = "#' + config.project_details['colour_hex_primary'] + '";</script>'
            extra_head += '<script>var secondary_colour = "#' + config.project_details['colour_hex_secondary'] + '";</script>'
            extra_head += '<script type="text/javascript" src="../keyword_history.js"></script>'
            extra_head += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

            html = ch.build_common_head("../../", extra_head)
            html += ch.build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; Tags &gt; ' + this_tag + '</p>', "../../")

            html += '<h1 id="pagetitle">Tag - ' + this_tag + '</h1>'
            html += '<h2>Tag History</h2>'

            # Place holders for the charts
            html += '<div id="papers_chart_div"></div>'
            html += '<div id="citations_chart_div"></div>'

            # List publications
            html += '<h2>Publications</h2>'
            html += '<p><em>' + str(len(zotero_tags[this_tag])) + ' publications with this tag</em></p>'

            fo.write(html)


            # Build the text needed for each paper
            #for this_paper in zotero_tags[this_tag]:
            for this_paper in utils.sort_hashes_by(papers, zotero_tags[this_tag], 'year'):

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

            html = ch.build_common_foot("../../")
            fo.write(html)

        fo.close()

    ######################################
    # Make HTML index page for all tags

    html = ch.build_common_head("../", "")
    html += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Tags</p>', "../")

    html += '<h1 id="pagetitle">Tags</h1>'
    html += '<p>' + str(len(zotero_tags)) + ' tags</p>'

    html_file.write(html)

    # Make a page with the tags on it with the number of papers and citations
    html_file.write('<table>')
    html_file.write('<tr><th style="text-align:left">Tag</th><th style="text-align:right">Number of papers</th><th style="text-align:right">Number of citations*</th></tr>')
    for this_tag in sorted(zotero_tags):

        # Some tags have a / in them, which messes with URL and file paths. Just swap it out for a - here.
        if this_tag.find("/") > 0:
            this_tag_safe = this_tag.replace("/", "-")
        else:
            this_tag_safe = this_tag


        #html = '<tr><td style="text-align:left"><a href="../tags/' + this_tag.replace(" ", "%20") + '/index.html">' + this_tag + '</a></td>'
        html = '<tr><td style="text-align:left"><a href="../tags/' + this_tag_safe + '/index.html">' + this_tag + '</a></td>'
        html += '<td style="text-align:right">' + str(zotero_tags_counts[this_tag]['total_papers']) + '</td>'
        html += '<td style="text-align:right">' + str(zotero_tags_counts[this_tag]['total_citations']) + '</td></tr>'
        html_file.write(html)
    html_file.write('</table>')

    html_file.write('<p>* Citation data from <a href="https://www.scopus.com">Scopus</a>.</p>')
    html_file.write('<p>Note that some publications may have multiple tags, this means that the total number of papers and citations may appear larger than for the project as a whole.</p>')

    html = ch.build_common_foot("../")
    html_file.write(html)
    ######################################
