
from config import config
from datetime import datetime
from . import utils
from . import common_html as ch

############################################################
# Home page with summary of years
############################################################
def build_home(papers):

    print("\n###HTML - Home###")

    summary = {}
    missing_year = {'num_papers': 0, 'citations': 0}

    html_file = open(config.html_dir + '/index.html', 'w', encoding='utf-8')
    data_file = open(config.html_dir + '/data.js', 'w', encoding='utf-8')

    temp = ch.build_common_head("", "")
    temp += ch.build_common_body("", "")
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
            except Exception:
                pass

        # If the year is not known, add it to the missing year list
        except Exception:
            missing_year['num_papers'] += 1
            try:
                missing_year['citations'] += int(this_paper['clean']['citations']['scopus']['count'])
            except Exception:
                pass

    # Add in some zeros when there are no papers for this year
    years = list(summary.keys())
    first_year = min(years)
    last_year = max(years)
    for this_year in range(int(first_year), int(last_year) + 1):
        try:
            summary[str(this_year)]['num_papers']
        except Exception:
            summary[str(this_year)] = {'num_papers': 0, 'cumulative': 0, 'citations': 0, 'cumulative_citations': 0}

    # Calculate the cumulative number of papers published
    for this_year in sorted(summary, reverse=False):
        try:
            summary[this_year]['cumulative'] = summary[this_year]['num_papers'] + summary[str(int(this_year) - 1)]['cumulative']
            summary[this_year]['cumulative_citations'] = summary[this_year]['citations'] + summary[str(int(this_year) - 1)]['cumulative_citations']
        except Exception:
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
    #cr_current_year = float(config.metrics_study_current_year)
    current_year = datetime.now().year
    cr_sum = 0.0
    cr_data_from = 0

    # Make a page with the headings on it
    html_file.write('<table>')
    html_file.write('<tr><th>Year</th><th>Number published</th><th>Cumulative</th><th>Citations* for papers published in this year</th><th>Cumulative citations* for papers published in this year</th></tr>')
    for this_year in sorted(summary, reverse=True):
        # Skip the years where nothing was published
        if summary[this_year]['num_papers'] == 0:
            continue

        cr_year = float(current_year - float(this_year))
        if cr_year < 1:
            cr_year = 1.0

        cr_sum += float(summary[this_year]['citations']) / cr_year
        cr_data_from += summary[this_year]['num_papers']

        # Build the table
        temp = '<tr><td><a href="papers/' + this_year + '/index.html">' + str(this_year) + '</a></td>'
        temp += '<td>' + utils.intWithCommas(summary[this_year]['num_papers']) + '</td>'
        temp += '<td>' + str(summary[this_year]['cumulative']) + '</td>'
        temp += '<td>' + utils.intWithCommas(summary[this_year]['citations']) + '</td>'
        temp += '<td>' + utils.intWithCommas(summary[this_year]['cumulative_citations']) + '</td></tr>'
        html_file.write(temp)

    # Unknown Row
    if missing_year['num_papers'] > 0:
        temp = '<tr>'
        temp += '<td style="font-size:12px;font-weight:bold;"><a href="papers/unknown/index.html">UNKNOWN</a></td>'
        temp += '<td>' + utils.intWithCommas(missing_year['num_papers']) + '</td>'
        temp += '<td>-</td>'
        temp += '<td>' + utils.intWithCommas(missing_year['citations']) + '</td>'
        temp += '<td>-</td>'
        temp += '</tr>'
        html_file.write(temp)
    html_file.write('</table>')

    temp = "<p>Publication year known for " + utils.intWithCommas(cr_data_from) + " of " + utils.intWithCommas(len(papers)) + " publications. <span class='help_text'>(<a href='about/index.html#missing_data'>What does this mean?</a>)</span></p>"
    temp += '<p>* Citation data from <a href="https://www.scopus.com">Scopus</a>.</p>'

    temp += ch.build_common_foot("./")
    html_file.write(temp)

    cr_sum = cr_sum / len(papers)
    return cr_sum, cr_data_from

