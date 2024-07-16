from config import config

###########################################################
# CSS colour scheme
###########################################################
def build_css_colour_scheme():

    # This function generates the CSS used for the colour scheme for the whole site.
    # This is includes the images for the top bar and naviagation bar.
    # The colours and image data is taken from the config file.

    html_file = open(config.html_dir + '/css/colour_scheme.css', 'w', encoding='utf-8')

    temp = "#header-container {background: #" + config.project_details['colour_hex_primary'] + ";}"

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

    temp += "#header-container:before,"
    temp += "#header-container:after {"
    temp += "display: block;"
    temp += "visibility: visible;"
    temp += "}"

    temp += ".nav-users a {"
    temp += "color: #fff;"
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
