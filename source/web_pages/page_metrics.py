import shutil

from config import config
from . import utils
from . import common_html as ch

###########################################################
# Build metrics page
###########################################################
def build_metrics(papers, age_weighted_citation, age_weighted_citation_data):

    print("\n###HTML - Metrics###")

    html_file = open(config.html_dir + '/metrics/index.html', 'w', encoding='utf-8')

    # Metric calculations
    total_publications = len(papers)
    total_citations = 0
    total_citations_data_from_count = 0
    paper_citations = []
    c20_index = 0
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
        except Exception:
            pass

    # We might not have any citations - e.g. if quotas hit.
    if total_citations > 0 and total_citations_data_from_count > 0:
        average_citations = float(total_citations)/float(total_citations_data_from_count)
    else:
        average_citations = 0

    # calculate h-index
    paper_citations.sort(reverse=True)
    h_index = 0

    for x in range(0, len(paper_citations)):
        if x > paper_citations[x]:
            break
        h_index = x

    # NUMBER OF PAPERS PER CITATION COUNT
    citation_number_limit = 100
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
        except Exception:
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
        except Exception:
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
        #if this_n_citations == round(average_citations, 0):
        #    colour = "#" + config.project_details['colour_hex_secondary']
        #if this_n_citations == median_citations:
        #    colour = "green"

        try:
            n_papers_with_x_citations += ",[" + str(this_n_citations) + "," + str(num_papers_citations[this_n_citations]) + ",'" + colour + "']"
        except Exception:
            n_papers_with_x_citations += ",[" + str(this_n_citations) + ",0,'" + colour + "']"

    n_papers_with_x_citations += "]);"

    # High Citations Range, but only if they are above the citation_number_limit
    plot_high_citation_chart = False
    if max_citations >= citation_number_limit:
        plot_high_citation_chart = True
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
                except Exception:
                    pass

            n_papers_with_x_citations += ",['" + str(bin_start) + "-" + str(bin_end) + "'," + str(num_papers_in_bin) + ",'" + colour + "']"

        n_papers_with_x_citations += "]);"

    # Put html together for this page

    shutil.copyfile(config.template_dir + '/metrics.js', config.html_dir + '/metrics/metrics.js')

    extra_head = '<script type="text/javascript" src="https://www.google.com/jsapi"></script>'
    extra_head += '<script type="text/javascript" src="../data.js"></script>'
    extra_head += '<script>' + n_papers_with_x_citations + '</script>'
    extra_head += '<script>var plot_high_citation_chart = "' + str(plot_high_citation_chart) + '";</script>'
    extra_head += '<script>var primary_colour = "#' + config.project_details['colour_hex_primary'] + '";</script>'

    extra_head += '<script type="text/javascript" src="metrics.js"></script>'

    temp = ch.build_common_head("../", extra_head)
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Metrics</p>', "../")

    temp += '<h1 id="pagetitle">Metrics</h1>'

    temp += '(All citation calculations based on <a href="https://www.scopus.com">Scopus</a> data.)'

    # Output Metrics
    temp += "<div class='metric_con'>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Total Publications</div>"
    temp += "<div class='metric_value'>" + utils.intWithCommas(total_publications) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(total_publications) + " Publications</div>"
    temp += "<div class='metric_description'>This is the number of all publications for the study.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Total Citations</div>"
    temp += "<div class='metric_value'>" + utils.intWithCommas(total_citations) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>This is the number of citations to all publications for the project.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Mean Citations Per Publication</div>"
    temp += "<div class='metric_value'>" + str("{0:.1f}".format(round(average_citations, 2))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The total number of citations divided by the total number of publications.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Age-weighted Mean Citations Per Publication</div>"
    temp += "<div class='metric_value'>" + str("{0:.1f}".format(round(age_weighted_citation, 3))) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(age_weighted_citation_data) + " Publications</div>"
    temp += "<div class='metric_description'>Age-weighted Mean Citations Per Publication.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>Median Citations</div>"
    temp += "<div class='metric_value'>" + str(median_citations) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The median number of citations for the project.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>h-index</div>"
    temp += "<div class='metric_value'>" + str(h_index) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>h-index is the largest number h such that h publications from a study have at least h citations.</div>"
    temp += "</div>"

    temp += "<div class='metric'>"
    temp += "<div class='metric_name'>c" + str(c_index_bound) + "-index</div>"
    temp += "<div class='metric_value'>" + utils.intWithCommas(c20_index) + "</div>"
    temp += "<div class='metric_stats_data'>Data From " + utils.intWithCommas(total_citations_data_from_count) + " Publications</div>"
    temp += "<div class='metric_description'>The number of publications that have at least " + str(c_index_bound) + " citations.</div>"
    temp += "</div>"

    #temp += "<div class='metric'>"
    #temp += "</div>"

    temp += "</div>"
    temp += "<div class='clear'></div>"

    temp += '<div id="cumulative_div"></div>'
    temp += "<p style='text-align:center;'>Data from " + utils.intWithCommas(age_weighted_citation_data) + " publications. <span class='help_text'>(<a href='../about/index.html#missing_data'>What does this mean?</a>)</span></p>"

    temp += '<div id="papers_per_year_div"></div>'
    temp += "<p style='text-align:center;'>Data from " + utils.intWithCommas(age_weighted_citation_data) + " publications. <span class='help_text'>(<a href='../about/index.html#missing_data'>What does this mean?</a>)</span></p>"

    temp += '<div id="papers_per_citation_count_div"></div>'
    #temp += "<div style='margin-left:auto;margin-right:auto;'><div class='average_citations' style='height:15px; width:33px; float:left; background:#" + config.project_details['colour_hex_secondary'] + "'></div><div style='height: 15px;line-height: 15px;padding-left: 40px;'> Mean number of citations</div></div>"
    #temp += "<div style='margin-left:auto;margin-right:auto;margin-top:5px;'><div class='average_citations' style='height:15px; width:33px; float:left; background:green'></div><div style='height: 15px;line-height: 15px;padding-left: 40px;'> Median number of citations</div></div>"
    temp += "<p style='text-align:center;'>Data from " + utils.intWithCommas(total_citations_data_from_count) + " publications. <span class='help_text'>(<a href='../about/index.html#missing_data'>What does this mean?</a>)</span></p>"

    if plot_high_citation_chart:
        temp += '<div id="papers_per_high_citation_count_div"></div>'
        temp += "<p style='text-align:center;'>Data from " + utils.intWithCommas(total_citations_data_from_count) + " publications. <span class='help_text'>(<a href='../about/index.html#missing_data'>What does this mean?</a>)</span></p>"

    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)
