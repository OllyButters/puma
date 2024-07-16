import os
import sys

from config import config
from . import common_html as ch

SITE_SECOND_TITLE = " Publications"

############################################################
# Home page with summary of years
############################################################
def build_years(papers):

    print("\n###HTML papers list###")

    yearly_papers = {}
    html_file = open(config.html_dir + '/papers/index.html', 'w', encoding='utf-8')

    temp = ch.build_common_head("../", "")
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Papers List</p>', "../")

    temp += '<h1 id="pagetitle">Papers List</h1>'
    # Altmetric include
    temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

    html_file.write(temp)
    main = "<p>"

    # Build the text needed for each paper
    for this_paper in papers:
        try:
            # Call draw paper function
            html = ch.draw_paper(this_paper, "../../")

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
        year_file = open(config.html_dir + '/papers/' + this_year + '/index.html', 'w', encoding='utf-8')

        temp = ch.build_common_head("../../", "")
        temp += ch.build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; <a href="../index.html">Papers List</a> &gt; ' + this_year + '</p>', "../../")

        temp += '<h1 id="pagetitle">Papers List - ' + this_year + '</h1>'
        temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

        temp += '<h2>' + str(len(yearly_papers[this_year])) + ' Publications From ' + this_year + '</h2>'
        year_file.write(temp)
        # This is a list
        for this_item in yearly_papers[this_year]:
            temp = list(this_item.values())
            year_file.write(temp[0])

        temp = ch.build_common_foot("../../")
        year_file.write(temp)

    main += "</p>" + ch.build_common_foot("../../")
    html_file.write(main)

    # == Publications from unknown years ==
    if not os.path.exists(config.html_dir + '/papers/unknown'):
        os.mkdir(config.html_dir + '/papers/unknown')
    unknown_file = open(config.html_dir + '/papers/unknown/index.html', 'w', encoding='utf-8')

    # Put html together for this page
    temp = '<!DOCTYPE html><html lang="en-GB">'

    # html head
    temp += '<head>'
    temp += '<title>' + SITE_SECOND_TITLE + '</title>'
    temp += '<link rel="stylesheet" href="../../css/style_main.css">'
    temp += '<link rel="stylesheet" href="../../css/colour_scheme.css">'
    temp += '</head>'

    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../../index.html">Home</a> &gt; <a href="../index.html">Papers List</a> &gt; Unknown Year</p>', "../../")

    temp += '<h1 id="pagetitle">Papers List - Unknown Year</h1>'
    # Altmetric include
    temp += "<script type='text/javascript' src='https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js'></script>"

    n = 0
    html = ""
    for this_paper in papers:
        try:
            this_paper['clean']['clean_date']['year']
        except:
            html += ch.draw_paper(this_paper, "../../")
            n += 1

    temp += '<h2>' + str(n) + ' Publications From Unknown Years</h2>'
    temp += '<h3>Some data may be missing from these publications.</h3>'
    temp += html

    temp += ch.build_common_foot("../")

    unknown_file.write(temp)

