import config.config as config
import shutil
import math


################################################################################
# Build an HTML report of the status of the important fields. Useful for cleaning
################################################################################
def coverage_report(papers):

    print('\n### Building coverage report ###')

    cov_css = '''
                .missing_required {background-color: red;}
                .missing_good_to_have {background-color: orange;}
                .raw_missing {color: red;}
                .raw_present {color: green;}
                tr:nth-child(even) {background-color: lightgray;}
                th {background-color: #4CAF50; color: white;}
                td, th {padding: 0.2em;}
                a:visited {color: red;}
    '''

    cov_html = '<table class="tablesorter">'
    cov_html += '''<thead><tr>
                    <th></th>
                    <th>Data</th>
                    <th>Zotero &uarr;&darr;</th>
                    <th>DOI &uarr;&darr;</th>
                    <th>PMID &uarr;&darr;</th>
                    <th>Scopus<br/>ID &uarr;&darr;</th>
                    <th>PMID DOI<br/>Lookup</th>
                    <th>PMID title<br/>Lookup</th>
                    <th>Title<br/>&uarr;&darr;</th>
                    <th>Keywords<br/>MeSH<br/>&uarr;&darr;</th>
                    <th>Keywords<br/>Other<br/>&uarr;&darr;</th>
                    <th>Abstract<br/>&uarr;&darr;</th>
                    <th>First<br/>Author &uarr;&darr;</th>
                    <th>Raw<br/>affil</th>
                    <th>First<br/>Author<br/>affil &uarr;&darr;</th>
                    <th>Clean<br/>Inst.<br/>&uarr;&darr;</th>
                    <th>Geocoded<br/>&uarr;&darr;</th>
                    <th>Country<br/>&uarr;&darr;</th>
                    <th>Clean<br/>Date &uarr;&darr;</th>
                    <th>Journal<br/>&uarr;&darr;</th>
                    <th>Volume<br/>&uarr;&darr;</th>
                    <th>Issue<br/>&uarr;&darr;</th>
                    <th>Scopus<br/>Citations &uarr;&darr;</th>
                    <th>Scopus<br/>Type &uarr;&darr;</th>
                </tr></thead>
                <tbody>'''

    status = {}
    status['hash'] = 0
    status['zotero'] = 0
    status['doi'] = 0
    status['pmid'] = 0
    status['scopus'] = 0
    status['title'] = 0
    status['keywords_mesh'] = 0
    status['keywords_other'] = 0
    status['abstract'] = 0
    status['first_author'] = 0
    status['first_author_affiliation'] = 0
    status['clean_institution'] = 0
    status['geocoded'] = 0
    status['country'] = 0
    status['clean_date'] = 0
    status['journal'] = 0
    status['volume'] = 0
    status['issue'] = 0
    status['scopus_citation'] = 0

    for this_paper in papers:
        cov_html += '\n<tr class="item">'

        # Put a non-functional checkbox in the first column to help with selecting papers
        cov_html += '<td><input type="checkbox"></td>'

        #####
        # Filename hash - this has to be present!
        try:
            fn_hash = this_paper['IDs']['zotero'] + '.cleaned.json'
            status['hash'] = status['hash'] + 1
        except:
            fn_hash = '???'

        cov_html += '<td><a href="processed/cleaned/' + fn_hash + '" target="_blank">data</a></td>'

        #####
        # Zotero ID - this has to be present!
        try:
            zotero = this_paper['IDs']['zotero']
            cov_html += '<td><a href = "http://www.zotero.org/groups/' + config.zotero_id + '/items/itemKey/' + zotero + '" target="_blank">' + zotero + '</a></td>'
            status['zotero'] = status['zotero'] + 1
        except:
            cov_html += '<td class="missing_required">???</td>'

        #####
        # DOI - Not required, but REALLY useful
        try:
            doi = this_paper['IDs']['DOI']
            if doi != '':
                cov_html += '<td><a href="https://doi.org/' + doi + '" target="_blank">' + doi + '</a></td>'
                status['doi'] = status['doi'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # PMID - Not required, but REALLY useful
        try:
            pmid = this_paper['IDs']['PMID']
            if pmid != '':
                cov_html += '<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/' + pmid + '" target="_blank">' + pmid + '</a></td>'
                status['pmid'] = status['pmid'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # scopus - Not required, but REALLY useful
        try:
            scopus_id = this_paper['IDs']['scopus']
            if scopus_id != '':
                cov_html += '<td><a href="https://api.elsevier.com/content/article/eid/' + scopus_id + '" target="_blank">' + scopus_id + '</a></td>'
                status['scopus'] = status['scopus'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # PMID lookup based on DOI
        try:
            doi = this_paper['IDs']['DOI']
            if doi != '':
                # Swap the / for the html encoded version
                doi.replace('/', '%2F')
                cov_html += '<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=' + doi + '" target="_blank">Do lookup</a></td>'
            else:
                raise Exception()
        except:
            # Not bothered if not there
            cov_html += '<td></td>'

        #####
        # PMID lookup based on title
        try:
            title = this_paper['clean']['title']
            if title != '':
                ascii_title = this_paper['clean']['title']
                cov_html += '<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=' + ascii_title + '" target="_blank">Do lookup</a></td>'
            else:
                raise Exception()
        except:
            # Not bothered if not there
            cov_html += '<td></td>'

        #####
        # Paper title
        try:
            title = this_paper['clean']['title']
            if title != '':
                cov_html += '<td title="' + title + '">OK</td>'
                status['title'] = status['title'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Keywords MeSH
        try:
            keywords_mesh = this_paper['clean']['keywords']['mesh']
            if keywords_mesh != '':
                if len(keywords_mesh) > 0:
                    # cov_html += '<td>OK</td>'
                    cov_html += '<td>' + str(len(keywords_mesh)) + '</td>'
                    status['keywords_mesh'] = status['keywords_mesh'] + 1
                else:
                    raise Exception()
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Keywords other
        try:
            keywords_other = this_paper['clean']['keywords']['other']
            if keywords_other != '':
                if len(keywords_other) > 0:
                    # cov_html += '<td>OK</td>'
                    cov_html += '<td>' + str(len(keywords_other)) + '</td>'
                    status['keywords_other'] = status['keywords_other'] + 1
                else:
                    raise Exception()
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Abstract
        try:
            abstract = this_paper['clean']['abstract']
            if abstract != '':
                cov_html += '<td>OK</td>'
                status['abstract'] = status['abstract'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # First author - Not required, but REALLY useful
        try:
            first_author = this_paper['clean']['first_author']
            if first_author != '':
                cov_html += '<td title = "' + first_author + '">OK</td>'
                status['first_author'] = status['first_author'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # RAW first author affiliation coverage
        try:
            raw_affil_pmid = this_paper['raw']['pmid_data']['MedlineCitation']['Article']['AuthorList'][0]['AffiliationInfo'][0]['Affiliation']
            raw_affil_pmid = '<span class="raw_present">P</span>'
        except:
            raw_affil_pmid = '<span class="raw_missing">P</span>'

        try:
            raw_affil_doi = this_paper['raw']['doi_data']['author'][0]['affiliation'][0]['name']
            raw_affil_doi = '<span class="raw_present">D</span>'
        except:
            raw_affil_doi = '<span class="raw_missing">D</span>'

        try:
            raw_affil_scopus = this_paper['raw']['scopus_data']['search-results']['entry'][0]['affiliation'][0]['affilname']
            raw_affil_scopus = '<span class="raw_present">S</span>'
        except:
            raw_affil_scopus = '<span class="raw_missing">S</span>'

        cov_html += '<td>' + raw_affil_pmid + ' ' + raw_affil_doi + ' ' + raw_affil_scopus + '</td>'

        #####
        # First author affiliation - Not required, but REALLY useful
        try:
            first_author_affiliation = this_paper['clean']['location']['candidate_institute']
            if first_author_affiliation != '':
                cov_html += '<td title = "' + first_author_affiliation + '">OK</td>'
                status['first_author_affiliation'] = status['first_author_affiliation'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # CLEAN first author affiliation - Not required, but REALLY useful
        try:
            clean_institution = this_paper['clean']['location']['clean_institute']
            if clean_institution != '':
                cov_html += '<td title = "' + clean_institution + '">OK</td>'
                status['clean_institution'] = status['clean_institution'] + 1
            else:
                raise Exception()
        except:
            # Check to see if we have a DIRTY first author affiliation, if we do then we really should have a clean one too.
            try:
                first_author_affiliation = this_paper['clean']['location']['candidate_institute']
                if first_author_affiliation != '':
                    cov_html += '<td class="missing_required">???</td>'
                else:
                    raise Exception()
            except:
                cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Geocoded
        try:
            latitude = this_paper['clean']['location']['latitude']
            longitude = this_paper['clean']['location']['longitude']
            if latitude != '' and longitude != '':
                cov_html += '<td>OK</td>'
                status['geocoded'] = status['geocoded'] + 1
            else:
                raise Exception()
        except:
            # Check to see if we have a clean institute, if we do then we should have a geocode
            try:
                clean_institution = this_paper['clean']['location']['clean_institute']
                if clean_institution != '':
                    cov_html += '<td class="missing_required">???</td>'
                else:
                    raise Exception()
            except:
                cov_html += '<td class="missing_good_to_have">???</!d>'

        #####
        # Country
        try:
            country = this_paper['clean']['location']['country']
            if country != '':
                cov_html += '<td>' + country + '</td>'
                status['country'] = status['country'] + 1
            else:
                raise Exception()
        except:
            # Check to see if we have a clean institute, if we do then we should have a geocode
            try:
                clean_institution = this_paper['clean']['location']['clean_institute']
                if clean_institution != '':
                    cov_html += '<td class="missing_required">???</td>'
                else:
                    raise Exception()
            except:
                cov_html += '<td class="missing_good_to_have">???</!d>'

        #####
        # Clean date
        try:
            clean_date = this_paper['clean']['clean_date']['year']
            if clean_date != '':
                cov_html += '<td>' + this_paper['clean']['clean_date']['year'] + '</td>'
                status['clean_date'] = status['clean_date'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_required">???</!d>'

        #####
        # Journal
        try:
            journal = this_paper['clean']['journal']['journal_name']
            if journal != '':
                cov_html += '<td>OK</td>'
                status['journal'] = status['journal'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Volume
        try:
            volume = this_paper['clean']['journal']['volume']
            if volume != '':
                cov_html += '<td>OK</td>'
                status['volume'] = status['volume'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Issue
        try:
            issue = this_paper['clean']['journal']['issue']
            if issue != '':
                cov_html += '<td>OK</td>'
                status['issue'] = status['issue'] + 1
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        #####
        # Scopus citations
        try:
            scopus_url = 'http://api.elsevier.com/content/search/scopus'

            # Start a blank string to add to so can gracefully bail if problems
            scopus_html = ''

            # These should exist, they might be empty though
            doi = this_paper['IDs']['DOI']
            pmid = this_paper['IDs']['PMID']

            try:
                # Might not exist
                scopus = str(this_paper['clean']['citations']['scopus']['count'])

                # If there is some scopus data then show it
                if scopus != '':
                    scopus_html += '<td>'
                    scopus_html += scopus + '&nbsp;'
                    status['scopus_citation'] = status['scopus_citation'] + 1
            except:
                scopus_html += '<td class="missing_good_to_have">'

            # Put links in
            if doi != '':
                scopus_html += '<a href="' + scopus_url + '?apiKey=' + config.scopus_api_key + '&query=DOI(' + doi + ')&httpAccept=application%2Fjson" target="_blank">DOI</a>'
            if pmid != '':
                # Stick in a space if needed
                if doi != '':
                    scopus_html += '&nbsp;'
                scopus_html += '<a href="' + scopus_url + '?apiKey=' + config.scopus_api_key + '&query=PMID(' + pmid + ')&httpAccept=application%2Fjson" target="_blank">PMID</a>'

            scopus_html += '</td>'
            cov_html += scopus_html
        except:
            cov_html += '<td class="missing_required">ERROR</td>'

        #####
        # Scopus type
        try:
            scopus_type = this_paper['raw']['scopus_data']['search-results']['entry'][0]['subtypeDescription']
            if scopus_type != '':
                cov_html += '<td>' + scopus_type + '</td>'
            else:
                raise Exception()
        except:
            cov_html += '<td class="missing_good_to_have">???</td>'

        cov_html += '</tr>'
    cov_html += '</tbody></table>'

    #####
    # Build a status table
    status_table = '<table>'
    status_table += '''<tr>
                    <th></th>
                    <th>Hash</th>
                    <th>Zotero</th>
                    <th>DOI</th>
                    <th>PMID</th>
                    <th>Scopus</th>
                    <th>Title</th>
                    <th>Keywords<br/>MeSH</th>
                    <th>Keywords<br/>Other</th>
                    <th>Abstract</th>
                    <th>First<br/>Author</th>
                    <th>First<br/>Author<br/>affil</th>
                    <th>Clean<br/>Inst</th>
                    <th>Geocoded</th>
                    <th>Country</th>
                    <th>Clean<br/>Date</th>
                    <th>Journal</th>
                    <th>Volume</th>
                    <th>Issue</th>
                    <th>Scopus<br/>Citations</th>
                </tr>'''

    number_of_papers = len(papers)

    # Actual values
    status_table += '<tr>'
    status_table += '<td>Number (out of ' + str(number_of_papers) + ')</td>'
    status_table += '<td>' + str(status['hash']) + '</td>'
    status_table += '<td>' + str(status['zotero']) + '</td>'
    status_table += '<td>' + str(status['doi']) + '</td>'
    status_table += '<td>' + str(status['pmid']) + '</td>'
    status_table += '<td>' + str(status['scopus']) + '</td>'
    status_table += '<td>' + str(status['title']) + '</td>'
    status_table += '<td>' + str(status['keywords_mesh']) + '</td>'
    status_table += '<td>' + str(status['keywords_other']) + '</td>'
    status_table += '<td>' + str(status['abstract']) + '</td>'
    status_table += '<td>' + str(status['first_author']) + '</td>'
    status_table += '<td>' + str(status['first_author_affiliation']) + '</td>'
    status_table += '<td>' + str(status['clean_institution']) + '</td>'
    status_table += '<td>' + str(status['geocoded']) + '</td>'
    status_table += '<td>' + str(status['country']) + '</td>'
    status_table += '<td>' + str(status['clean_date']) + '</td>'
    status_table += '<td>' + str(status['journal']) + '</td>'
    status_table += '<td>' + str(status['volume']) + '</td>'
    status_table += '<td>' + str(status['issue']) + '</td>'
    status_table += '<td>' + str(status['scopus_citation']) + '</td>'
    status_table += '</tr>'

    # Percentages
    status_table += '<tr>'
    status_table += '<td>Percentage</td>'
    status_table += '<td>' + str(int(math.floor(100*status['hash']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['zotero']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['doi']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['pmid']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['scopus']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['title']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['keywords_mesh']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['keywords_other']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['abstract']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['first_author']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['first_author_affiliation']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['clean_institution']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['geocoded']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['country']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['clean_date']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['journal']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['volume']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['issue']/number_of_papers))) + '</td>'
    status_table += '<td>' + str(int(math.floor(100*status['scopus_citation']/number_of_papers))) + '</td>'
    status_table += '</tr></table>'

    # Title
    title = '<h1>' + config.project_details['name'] + ' coverage report</h1>'

    # Key
    key = '<table><tr><td class="missing_good_to_have">&nbsp;</td><td>Missing</td></tr><tr><td class="missing_required">&nbsp;</td><td>Required</td></tr></table>'

    scripts = '<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>'
    scripts += '<script type="text/javascript" src="jquery.tablesorter.js"></script>'
    scripts += "<script>$(function(){$('table').tablesorter({widgets        : ['zebra', 'columns'],usNumberFormat : false,sortReset      : true,sortRestart    : true});});</script>"

    output_text = '<html><head><style>' + cov_css + '</style>' + scripts + '</head><body>' + title + status_table + key + '\n<br/><br/>' + cov_html + '</body></html>'
    coverage_file = open(config.cache_dir + '/coverage_report.html', 'w')
    coverage_file.write(output_text)

    # Copy the jquery.tablesorter.js file across
    shutil.copy(config.template_dir + '/jquery.tablesorter.js', config.cache_dir + '/jquery.tablesorter.js')
