from config import config
from . import common_html as ch

###########################################################
# Help Page
###########################################################
def build_about():

    print("\n###HTML - about page###")

    html_file = open(config.html_dir + '/about/index.html', 'w', encoding='utf-8')

    # # Put html together for this page

    temp = ch.build_common_head("../", "")
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; About</p>', "../")

    temp += '<h1 id="pagetitle">About</h1>'

    temp += '<ul>'
    temp += '<li><a href="#where_does_the_data_come_from">Where does the data come from?</a></li>'
    temp += '<li><a href="#what_des_pubmed_doi_and_api_mean">What does PubMed, DOI and API mean?</a></li>'
    temp += '<li><a href="#where_does_citation_data_come_from">Where does the citation data come from?</a></li>'
    temp += '<li><a href="#geo_data">Geolocation data</a></li>'
    temp += '<li><a href="#missing_data">Why don\'t some statistics use data from all publications?</a></li>'
    temp += '<li><a href="#altmetric_circles">What are the circles with numbers in them?</a></li>'
    temp += '<li><a href="#all_keywords_vs_mesh">What\'s the difference between "All keywords" and "Major keywords (MeSH)"?</a>>/li>'
    temp += '<li><a href="#missing_from_keywords">Why doesn\'t my paper show up when I click on a keyword which describes it?</a></li>'
    temp += '<li><a href="#publication_dates">What are some publications in a different year to what I would expect?</a></li>'
    temp += '<li><a href="#what_is_a_word_cloud">What is a word cloud?</a></li>'
    temp += '<li><a href="#i_want_to_know_more">I want to know more!</a></li>'
    temp += '<li><a href="#how_was_this_funded">How was this project funded?</a></li>'
    temp += '</ul>'

    temp += '<h2 id="where_does_the_data_come_from">Where does the data come from?</h2>'
    temp += '<p>It is mostly <em>metadata</em> we use here. You can think of the journal articles themselves as being the <em>data</em> part, with information about the journal articles (e.g. author lists, keywords, citations etc) being the <em>metadata</em> part. We get this metadata for the journal articles from the PubMed and DOI APIs.</p>'

    temp += '<h2 id="what_des_pubmed_doi_and_api_mean">What does PubMed, DOI and API mean?</h2>'
    temp += '<p>PubMed is an online search engine which contains journal articles on lifesciences and biomedical topics. DOI is short for Digital Object Identifier and is the de facto way of referring to journal articles. An example DOIs is <a href="https://doi.org/10.12688/f1000research.25484.2">https://doi.org/10.12688/f1000research.25484.2</a> and <a href="http://doi.org">http://doi.org</a> are in charge of maintaining the official list of them. An API is an Application Programming Interface, which is a fancy way of saying that computer programs can send and receive data to a service. So here we send and receive metadata to/from PubMed and doi.org.</p>'

    temp += '<h2 id="where_does_citation_data_come_from">Where does the citation data come from?</h2>'
    temp += '<p>Citation data is retrieved from the <a href="https://www.elsevier.com/solutions/scopus">Scopus</a> API provided by Elsevier. This is cached locally and updated regularly. You can see the last time the citation data was updated in the last updated note at the bottom of the page.</p>'

    temp += '<h2 id="geo_data">Geolocation data</h2>'
    temp += '<p>We use the first author\'s postal address and email address to assign a journal article to a location. Sometimes this may not exist in the metadata, which is why not all journal articles will be plotted on the maps. Occasionally a journal article indicates the authors wish a different author to be considered the lead author, we cannot process this information, so the first author is used.</p>'
    temp += '<p>The coordinate location of institutions is retrieved from the free open data project <a href="https://www.wikidata.org">Wikidata</a>.</p>'

    temp += '<h2 id="missing_data">Why don\'t some statistics use data from all publications?</h2>'
    temp += '<p>Metadata can be missing because some journal articles are very old and the metadata about it just doesn\'t exist, or sometimes a recent journal article may not have the metadata about it available <em>yet</em>. There are many different ways to track journal articles'
    temp += ' and prior to the introduction of DOIs in 2000 there was no standard method. This means some metadata on old journal articles could have been lost or not recorded.</p>'

    temp += '<p>The metadata used for the statistics is gathered from databases which only collect data from particular journals (no one has 100% of journal articles in their database). This means some statistics shown here are under-reported (e.g. citation counts and their derivatives) or may be missing entirely (e.g. journal articles with no location assigned to them).'

    temp += '<h2 id="altmetric_circles">What are the circles with numbers in them?</h2>'
    temp += '<p>Journal article citations are one way to track how an article is being used. Another way is by considering a wider set of metrics, such as if it is mentioned in news articles, or on social media. This is what <a href"https://www.altmetric.com">https://www.altmetric.com</a> does, and by hovering your mouse over one of the circles (or clicking on it) you can get an overview of where this article is being talked about.</p>'

    temp += '<h2 id="all_keywords_vs_mesh">What\'s the difference between "All keywords" and "Major keywords (MeSH)"?</h2>'
    temp += '<p>The <em>All keywords</em> page shows author defined keywords and MeSH terms. MeSH (Medical Subject Headings) is a commonly used hierarchical controlled vocabulary. This means that there is a list of well defined terms which are allowed to be used (controlled vocabulary), and they are related to one another (hierarchical) e.g. "finger" belongs to "hand". So the Major Keywords are the broader areas, and all keywords will be more fine grained. MeSH terms are applied retrospectively to a subset of all journal articles.</p>'

    temp += '<h2 id="missing_from_keywords">Why doesn\'t my paper show up when I click on a keyword which describes it?</h2>'
    temp += '<p>Authors can often add whatever they want for their keywords when they publish a journal article, this may result in slightly different words for the same thing being used by different articles. Sometimes the author defined keywords are not present in the metadata at all.</p>'

    temp += '<h2 id="publication_dates">What are some publications in a different year to what I would expect?</h2>'
    temp += '<p>How do <em>you</em> define when a journal article is published? Is it when it is accepted, or when a pre-print is available, or when it is made available on a journal\'s website, or when it appears in a journal\'s printed volume, or something else? Different journals have different preferences for which they prefer. This will likely mean that some journal articles appear later here than you may expect them to.</p>'

    temp += '<h2 id="what_is_a_word_cloud">What is a word cloud?</h2>'
    temp += '<p>A word cloud is a way of showing how often words appear in a list of words. The more frequent a word is the bigger the word appears.</p>'

    temp += '<h2 id="i_want_to_know_more">I want to know more!</h2>'
    temp += '<ul>'
    temp += '<li>The source code is at <a href="https://github.com/OllyButters/puma" target="_blank">https://github.com/OllyButters/puma</a>. You can download and run this yourself!</li>'
    temp += '<li>Some documentation is at <a href="https://github.com/OllyButters/puma/wiki" target="_blank">https://github.com/OllyButters/puma/wiki</a></li>'
    temp += '<li>There is even a paper written about it: <a href="https://f1000research.com/articles/9-1095/v2" target="_blank">https://f1000research.com/articles/9-1095/v2</a></li>'
    temp += '<li>You can talk to us at: <a href="https://twitter.com/DrOllyButters" target="_blank">@DrOllyButters</a>, <a href="https://twitter.com/DrBeccaWilson" target="_blank">@DrBeccaWilson</a> and <a href="https://twitter.com/_hugh_garner_" target="_blank">@_hugh_garner_</a></li></ul>'

    temp += '<h2 id="how_was_this_funded">How was this project funded?</h2>'
    temp += 'This project has been funded by:'
    temp += '<ul>'
    temp += '<li>CLOSER, whose mission is to maximise the use, value and impact of longitudinal studies. CLOSER is funded by the Economic and Social Research Council (ESRC) and Medical Research Council (MRC) (grant reference: ES/K000357/1).</li>'
    temp += '<li>Becca Wilson is a UKRI Innovation Fellow with HDR UK [MR/S003959/1].</li>'
    temp += '<li>The Nuffield Foundation research placement program.</li>'
    temp += '<li>The Wellcome Trust and Medical Research Council (grant number 108439/Z/15/Z).</li>'
    temp += '<li>The European Union’s Horizon 2020 research and innovation programme under grant agreement No 824989 and the Canadian Institutes of Health Research (CIHR).</li>'
    temp += '<li>The National Institute for Health Research Applied Research Collaboration.</li>'
    temp += '</ul>'

    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)
