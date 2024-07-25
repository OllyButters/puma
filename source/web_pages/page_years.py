import os

from config import config
from . import common_html as ch
from . import utils

SITE_SECOND_TITLE = " Publications"

############################################################
# Home page with summary of years
############################################################
def build_years(papers):

    print("\n###HTML papers list###")

    years = {}

    # Build a dict of years with a list of each hash (paper) that has this year.
    for this_paper in papers:
        try:
            # If this tag is not already in the dict then add it
            if this_paper['clean']['clean_date']['year'] not in years:
                years[this_paper['clean']['clean_date']['year']] = list()
            years[this_paper['clean']['clean_date']['year']].append(this_paper['IDs']['hash'])
        except Exception:
            if 'unknown' not in years:
                years['unknown'] = list()
            years['unknown'].append(this_paper['IDs']['hash'])


    # Output the info into an HTML file
    # For each year dict item
    for this_year in years.keys():
        # Check there is some data for this year - not all do
        if len(years[this_year]) == 0:
            continue

        if not os.path.exists(config.html_dir + '/papers/' + this_year):
            os.mkdir(config.html_dir + '/papers/' + this_year)
        year_file = open(config.html_dir + '/papers/' + this_year + '/index.html', 'w', encoding='utf-8')

        temp = ch.build_common_head("../../", "")
        temp += ch.build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; <a href="../index.html">Papers List</a> &gt; ' + this_year + '</p>', "../../")

        temp += '<h1 id="pagetitle">Papers List - ' + this_year + '</h1>'
        temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

        temp += '<h2>' + str(len(years[this_year])) + ' Publications From ' + this_year + '</h2>'
        year_file.write(temp)

        for this_paper in utils.sort_hashes_by(papers, years[this_year], 'year'):

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

                year_file.write(html)

            except Exception:
                pass

        temp = ch.build_common_foot("../../")
        year_file.write(temp)
