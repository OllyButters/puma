#!/usr/bin/env python3

import logging

import config.config as config


# Output a bibtex file of the publications
def bibtex(papers):
    print("### Building BibTeX file ###")

    articles = []

    for this_paper in papers:

        # read the cache
        try:
            # Entry type and name
            this_article = '@article{'
            this_article += this_paper['IDs']['hash']
            this_article += ',\n'

            # Author
            this_article += 'author = "'
            this_article += this_paper['clean']['full_author_list'][0]['family']
            this_article += '",\n'

            # Title
            this_article += 'title = "'
            this_article += this_paper['clean']['title']
            this_article += '",\n'

            # Journal
            this_article += 'journal = "'
            # this_artice['cleaned-journal'] set by clean.clean_journal
            this_article += this_paper['clean']['journal']['journal_name']
            this_article += '",\n'

            # Journal volume
            # this_article += 'volume = "'
            # this_article += papers[this_pmid]['JournalVolume']
            # this_article += '",\n'

            # Year
            this_article += 'year = "'
            this_article += this_paper['clean']['clean_date']['year']
            this_article += '",\n'

            # Close this one.
            this_article += '}\n'

            articles.append(this_article)

        except:
            logging.error("Cannot create bibtex output" + str(this_paper))
            pass

    # output to file in data_dir
    with open(config.data_dir + '/bibtex_list.bib', 'w') as bibfile:
        for a in articles:
            bibfile.write(a)
            bibfile.write('\n')
