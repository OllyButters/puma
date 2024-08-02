import os
import shutil
import time

from html import escape

from config import config

SITE_SECOND_TITLE = " Publications"


############################################################
# Draw Paper Function
# This function is used in the different paper lists to
# display a consistently formatted list.
############################################################
def draw_paper(this_paper, nav_path="./"):
    html = '<div class="paper">'

    # altmetric data
    try:
        if this_paper['IDs']['DOI']:
            html += '<div style="float:right; width:64px; height:64px;" data-badge-popover="left" data-badge-type="donut" data-doi="' + this_paper['IDs']['DOI'] + '" data-hide-no-mentions="true" class="altmetric-embed"></div>'
    except Exception:
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
        except Exception:
            pass
            #logging.debug('This author dropped from author list on webpage for some reason: ' + str(this_author))

    html += escape('; '.join(authors))
    html += '<br/>'

    try:
        if this_paper['clean']['journal']['journal_name'] != "":
            html += this_paper['clean']['journal']['journal_name']
    except Exception:
        pass

    try:
        if this_paper['clean']['journal']['volume'] != "":
            html += ', Volume ' + this_paper['clean']['journal']['volume']
    except Exception:
        pass

    try:
        if this_paper['clean']['journal']['issue'] != "":
            html += ', Issue ' + this_paper['clean']['journal']['issue']
    except Exception:
        pass

    try:
        if this_paper['clean']['clean_date']['year'] != "":
            html += " (" + str(this_paper['clean']['clean_date']['year']) + ")"
    except Exception:
        pass
    
    html += '<br/>'

    # PMID
    try:
        if this_paper['IDs']['PMID']:
            html += 'PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + str(this_paper['IDs']['PMID'])+'" target="_blank">' + str(this_paper['IDs']['PMID']) + '</a>&nbsp;'
    except Exception:
        pass

    # DOI
    try:
        if this_paper['IDs']['DOI']:
            html += 'DOI: <a href="https://doi.org/' + this_paper['IDs']['DOI'] + '" target="_blank">' + this_paper['IDs']['DOI'] + '</a>&nbsp;'
    except Exception:
        pass

    # Citation Counts and Sources
    html += '<br/>Citations: '

    try:
        # Try to display citation count with link to scopus page
        scopus_links = this_paper['raw']['scopus_data']['search-results']['entry'][0]['link']
        for this_link in scopus_links:
            if this_link['@ref'] == 'scopus-citedby':
                if this_paper['clean']['citations']['scopus']['count'] !='0':
                    html += '<a href="' + this_link['@href'] + '" target="_blank">' + str(this_paper['clean']['citations']['scopus']['count']) + '</a> (Scopus)'
                else:
                    html += '-'
    except Exception:
        try:
            html += str(this_paper['clean']['citations']['scopus']['count']) + ' (Scopus)'
        except Exception:
            html += '-'

    # Tags
    try:
        if len(this_paper['clean']['zotero_tags']) > 0:
            html += '<br/>Tags: '
            html_tags = []
            for this_tag in this_paper['clean']['zotero_tags']:
                html_tags.append('<a href="' + nav_path + 'tags/' + this_tag['tag'] + '/index.html">' + this_tag['tag'] + '</a>')
            html += '&nbsp;|&nbsp;'.join(html_tags)
    except Exception:
        pass

    html += '</div>'

    return html

############################################################
# Build common head for all pages
############################################################
def build_common_head(nav_path, extra_head_content):

    # Copy files across
    shutil.copyfile(config.template_dir + '/style_main.css', config.html_dir + '/css/style_main.css')
    shutil.copyfile(config.template_dir + '/timestamp.js', config.html_dir + '/timestamp.js')
    shutil.copyfile(config.template_dir + '/iframe.js', config.html_dir + '/iframe.js')

    # Put html together for this page
    head = '<!DOCTYPE html><html lang="en-GB">'

    head += '<head>'
    head += '<title>' + SITE_SECOND_TITLE + '</title>'
    head += '<meta charset="UTF-8">'
    head += '<link rel="stylesheet" href="' + nav_path + 'css/style_main.css">'
    head += '<link rel="stylesheet" href="' + nav_path + 'css/colour_scheme.css">'
    head += '<script src="' + nav_path + 'timestamp.js"></script>'
    head += '<script src="' + nav_path + 'iframe.js" defer></script>'

    head += extra_head_content
    head += '</head>'

    return head

############################################################
# Build common body for all pages
############################################################
def build_common_body(breadcrumb, nav_path):
    # This function builds the common header and nav bar for all pages.
    # nav_path used for changes to relative pathing depending on the page 
    # (i.e. Home does not need anything but next level down needs leading ../)

    html = '<body>'

    html += "<div id='header-container'>"
    html += "<div class='header width-master' role='banner'>"

    if os.path.isfile(config.config_dir + '/' + config.project_details['header_institution_logo_filename']):
        shutil.copy(config.config_dir + '/' + config.project_details['header_institution_logo_filename'], config.html_dir + '/' + config.project_details['header_institution_logo_filename'])
        html += '<div id="main-logo"><a accesskey="1" title="' + config.project_details['header_institution_name'] + '" href="' + config.project_details['header_institution_url'] + '"><img id="main-logo-img" src="' + nav_path + config.project_details['header_institution_logo_filename'] + '" alt="Institution logo"/></a></div>'

    html += "<div class='maintitle' id='maintitle1'>"
    html += "<span id='title1'><a href='" + nav_path + "index.html'>" + config.project_details['name'] + " - " + SITE_SECOND_TITLE + "</a></span>"
    html += "</div>"
    html += "</div>"
    html += "</div>"

    html += '<div id="wrapper" class="width-master">'
    html += '<div id="col1" role="navigation">'
    html += '<!-- navigation object : navigation -->'

    html += '<ul class="navgroup">'
    html += '<li><a href="' + nav_path + 'index.html">Home</a></li>'
    html += '<li><a href="' + nav_path + 'about/index.html">About</a></li>'
    if config.WEB_PAGE_PROJECT_ABOUT_HTML_FILE is not None:
        html += '<li><a href="' + nav_path + 'about/project.html">About project</a></li>'


    html += '<li><a href="' + nav_path + 'search/index.html">Search</a></li>'

    if config.WEB_PAGE_SHOW_ZOTERO_TAGS:
        html += '<li><a href="' + nav_path + 'tags/index.html">Tags</a></li>'

    html += '<li><a href="' + nav_path + 'keywords/index.html">All Keywords</a></li>'
    html += '<li><a href="' + nav_path + 'mesh/index.html">Major Keywords (MeSH)</a></li>'

    if config.WEB_PAGE_SHOW_INSTITUTE_COUNTRY_MAP:
        html += '<li><a href="' + nav_path + 'country/index.html">Map by Country</a></li>'

    if config.WEB_PAGE_SHOW_INSTITUTE_UK_MAP:
        html += '<li><a href="' + nav_path + 'institute/index.html">Map by UK institute</a></li>'

    html += '<li><a href="' + nav_path + 'metrics/index.html">Metrics</a></li>'
    html += '<li><a href="' + nav_path + 'keyword_wordcloud/index.html">Keyword Cloud</a></li>'
    html += '<li><a href="' + nav_path + 'abstractwordcloud/index.html">Abstract Word Cloud</a></li>'

    html += '</ul>'

    html += '<div class="after-navgroup">'
    html += '<!-- navigation object : navigation bottom -->'
    html += '<!-- start navigation : additional logo -->'

    html += '</div>'
    html += '</div>'

    html += breadcrumb

    html += '<div id="content" role="main">'

    return html


###########################################################
# Build a common footer for all the pages.
###########################################################
def build_common_foot(nav_path):

    # Output the timestamp to file
    timestamp_file = open(config.html_dir + '/timestamp.txt', 'w', encoding='utf-8')
    timestamp_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))
    timestamp_file.close()

    html = '</div>'
    html += '</div>'
    html += '<div class="foot">'
    # Put a link to the timestamp file, but update it to the content via JS. This means
    # there will be a mechanism to see this when running locally as COS stops JS calls.
    html += '<span id="update_timestamp"><a href="' + nav_path + '/timestamp.txt">Update time</a></span>'
    html += '<script type="text/javascript" src="' + nav_path + '/timestamp.js"></script>'
    html += '<script>update_timestamp("' + nav_path + '/timestamp.txt")</script>'
    html += '</div>'
    html += '</body>'
    html += '</html>'

    return html