
from config import config
from . import page_maps
from . import page_home
from . import page_wordclouds
from . import page_tags
from . import page_keywords
from . import page_about
from . import page_search
from . import page_css
from . import page_years
from . import page_metrics

def build_all(papers, papers_with_keywords, papers_with_abstract_text):
    age_weighted_citations, age_weighted_citations_data = page_home.build_home(papers)
    page_years.build_years(papers)
    page_maps.build_country_map(papers)
    page_maps.build_UK_institute_map(papers)
    page_wordclouds.build_keyword_word_cloud(papers_with_keywords)
    page_wordclouds.build_abstract_word_cloud(papers_with_abstract_text)
    page_keywords.build_mesh(papers)
    page_about.build_about()
    page_about.build_about_study()
    page_search.build_search(papers)
    page_css.build_css_colour_scheme()
    page_metrics.build_metrics(papers, age_weighted_citations, age_weighted_citations_data)

    if config.WEB_PAGE_SHOW_ZOTERO_TAGS:
        page_tags.build_zotero_tags(papers)
