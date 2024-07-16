import json
import shutil

from config import config
from . import common_html as ch


###########################################################
# Search Page
###########################################################
def build_search(papers):

    print("\n###HTML - Search page###")

    html_file = open(config.html_dir + '/search/index.html', 'w')
    shutil.copyfile(config.template_dir + '/search.js', config.html_dir + '/search/search.js')

    # Put html together for this page

    extra_head = '<script>var name = "' + config.project_details['name'] + '";</script>'
    extra_head += '<script src="search.js"></script>'
    extra_head += '<style> button { height:2em; font-size:1em; margin-left:5px; } input#search { width:50%; }</style>'

    temp = ch.build_common_head("../", extra_head)
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Search</p>', "../")

    temp += '<h1 id="pagetitle">Search</h1>'

    temp += '<p>Search for the fields titles, abstracts, keywords and authors.<br/>To narrow down the results search for multiple fields at once.</p>'
    temp += '<p><input type="text" id="search"><button onclick="search();">Search</button></p>'

    temp += '<div style="display:none;" id="searching">Searching...</div>'
    temp += '<h2 id="num_search_results"></h2>'
    temp += '<div id="search_results"></div>'

    searchable_data = []
    for this_paper in papers:
        this_subset = {}
        this_subset['IDs'] = {}
        
        try:
            this_subset['IDs']['DOI'] = this_paper['IDs']['DOI']
        except:
            pass
        
        try:
            this_subset['IDs']['PMID'] = this_paper['IDs']['PMID']
        except:
            pass
        
        #this_subset['clean'] = {}
        try:
            this_subset['title'] = this_paper['clean']['title']
        except:
            pass

        try:
            this_subset['abstract'] = this_paper['clean']['abstract']
        except:
            pass

        try:
            this_subset['keywords'] = this_paper['clean']['keywords']
        except:
            pass

        try:
            this_subset['full_author_list'] = []
            for this_author in this_paper['clean']['full_author_list']:
                this_subset['full_author_list'].append(this_author['clean'])
        except:
            pass

        try:
            this_subset['journal'] = this_paper['clean']['journal']
        except:
            pass

        searchable_data.append(this_subset)

    temp += '<script id="search_data">var papers = ' + str(json.dumps(searchable_data)).replace("<", "&lt;").replace(">", "&gt;") + ';</script>'

    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)
