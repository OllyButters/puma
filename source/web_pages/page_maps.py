import shutil

from config import config
from . import utils
from . import common_html as ch

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

    html_file = open(config.html_dir + '/country/index.html', 'w', encoding='utf-8')

    # Put html together for this page

    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/country/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/css/map.css')

    extra_head = '<link rel="stylesheet" href="../css/map.css">'
    extra_head += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script> <script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    extra_head += '<script type="text/javascript">' + "google.charts.load('current', {'packages':['geochart']});google.charts.setOnLoadCallback(drawRegionsMap);function drawRegionsMap() {var data = google.visualization.arrayToDataTable([ ['Country', 'Publications']" + country_string + "]); var options = { colorAxis: {colors: ['#" + config.project_details['colour_hex_secondary'] + "', '#" + config.project_details['colour_hex_primary'] + "']} }; var chart = new google.visualization.GeoChart(document.getElementById('regions_div')); chart.draw(data, options); }</script>"

    temp = ch.build_common_head("../", extra_head)
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications by Country</p>', "../")

    temp += '<h1 id="pagetitle">Publications by Country</h1>'

    temp += '<div id="regions_div" style="width: 900px; height: 500px;"><img src="loading.gif" alt="Loading"></div>'
    temp += "<p>Data from " + utils.intWithCommas(number_of_points) + " publications. <span class='help_text'>(<a href='about/index.html#missing_data'>What does this mean?</a>)</span></p>"

    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)


###########################################################
# Publications by UK institute
###########################################################
def build_UK_institute_map(papers):

    print("\n###HTML - Institute Map###")

    institutes = {}
    number_of_points = 0
    for this_paper in papers:
        try:
            if (
                this_paper['clean']['location']['clean_institute'] != "" and
                this_paper['clean']['location']['latitude'] != "" and 
                this_paper['clean']['location']['longitude'] != "" and
                    (
                    this_paper['clean']['location']['country'] == "United Kingdom" or
                    this_paper['clean']['location']['country'] == "Northern Ireland"
                    )
                ):
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

    html_file = open(config.html_dir + '/institute/index.html', 'w', encoding='utf-8')

    # # Put html together for this page

    shutil.copyfile(config.template_dir + '/loading.gif', config.html_dir + '/institute/loading.gif')
    shutil.copyfile(config.template_dir + '/map.css', config.html_dir + '/css/map.css')

    extra_head = '<link rel="stylesheet" href="../css/map.css">'
    extra_head += '<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>'
    extra_head += "<script>google.charts.load('current', {'packages':['geochart']});google.charts.setOnLoadCallback(drawMarkersMap);function drawMarkersMap() {var data = google.visualization.arrayToDataTable([['lat', 'lon', 'Institute','Publication count']" + institute_string + " ]); var options = {magnifyingGlass: {zoomFactor: '15.0'}, region: 'GB', displayMode: 'markers', colorAxis: {colors: ['#" + config.project_details['colour_hex_secondary'] + "', '#" + config.project_details['colour_hex_primary'] + "']}}; var chart = new google.visualization.GeoChart(document.getElementById('regions_div'));chart.draw(data, options); };</script>"

    temp = ch.build_common_head("../", extra_head)
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Publications by UK City</p>', "../")

    temp += '<h1 id="pagetitle">Publications by UK Institute</h1>'

    temp += '<div id="regions_div" style="width: 100%; min-width:500px; min-height: 500px;"><img src="loading.gif" alt="Loading"></div>'
    temp += "<p>Data from " + utils.intWithCommas(number_of_points) + " UK publications. <span class='help_text'>(<a href='../about/index.html#missing_data'>What does this mean?</a>)</span></p>"
    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)

