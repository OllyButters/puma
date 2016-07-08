#! /usr/bin/env python


def bibtex(papers):

    for this_paper in papers:

        # read the cache
        try:
            # Entry type and name
            this_article = '@article{'
            this_article += this_paper['IDs']['hash']
            this_article += ',\n'

            # Author
            this_article += 'author = "'
            this_article += this_paper['author'][0]['family']
            this_article += '",\n'

            # Title
            this_article += 'title = "'
            this_article += this_paper['title']
            this_article += '",\n'

            # Journal
            this_article += 'journal = "'
            this_article += this_paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
            this_article += '",\n'

            # Journal volume
            # this_article += 'volume = "'
            # this_article += papers[this_pmid]['JournalVolume']
            # this_article += '",\n'

            # Year
            this_article += 'year = "'
            this_article += this_paper['PubmedData']['History'][0]['Year']
            this_article += '",\n'

            # Close this one.
            this_article += '}\n'

        except:
            pass
