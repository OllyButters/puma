#! /usr/bin/env python

def bibtex(pmids, papers):
    
    import json
    import logging
        

    for this_pmid in pmids:

        #read the cache
        try:
            #Entry type and name
            this_article = '@article{'
            this_article += this_pmid
            this_article += ',\n'

            #Author
            this_article += 'author = "'
            this_article += papers[this_pmid]['AuthorList'][0]['LastName']
            this_article += '",\n'

            #Title
            this_article += 'title = "'
            this_article += papers[this_pmid]['ArticleTitle']
            this_article += '",\n'

            #Journal
            this_article += 'journal = "'
            this_article += papers[this_pmid]['Journal']
            this_article += '",\n'

            #Journal volume
            this_article += 'volume = "'
            this_article += papers[this_pmid]['JournalVolume']
            this_article += '",\n'

            #Year
            this_article += 'year = "'
            this_article += papers[this_pmid]['Year']
            this_article += '",\n'

            #Close this one.
            this_article += '}\n'

            #logging.info(str(this_pmid)+' not in citation cache')

            print this_article

        except:
            pass


