#Overview#

The general point of this 'pipeline' is to get meta-data about published papers and do some analysis/plotting of it. The general worflow is:

* Start with a list of papers (e.g. with a PubMed ID), and get the metadata for each.
* Build a cache of each bit of metadata.
* Pull out the required metadata from the cache into a cut down version.
* Clean it.
* Add some extra data to it (citations etc).
* Analyse it.
* Make some web pages showing it.
* Build a bibliographic reference of all of them for redistribution.


#Pre-requisites#

##Dependencies##

###Biopython###
To pull out any PubMed data the biopython library is needed:

sudo apt-get install python-biopython

main docs here
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch

http://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html

maybe useful info on field descriptions
http://www.nlm.nih.gov/bsd/mms/medlineelements.html

###Google docs###
Can interact with google docs - not currently using this.

sudo apt-get install python-gdata

##Source data##

Currently the whole process kicks off with a list of PubMed IDs. This is generated from a google docs spreadsheet and downloaded to the inputs directory.


#Program flow#

Run the program from /source/papers.py

##Get##


##Clean##
###Institute###
Names of institutes change all the time (especially over the lifetime of long studies like cohort studies). Then there are the many different spellings and formats of an institution name that people put on papers. In order to merge all of these together a CSV file is used that has all the different spellings mapped to a canonical institute name. It is a manual step to generate this list, and will probably require running the program a few times - adding to the list incrementally.

##Add##
###Geocoding###
Using the canonical institute name a geocode can be assigned - this is another CSV that has been manually generated.

##Analyse##

##HTML##


#Notes#
##Dates##
How do you define when a paper was published? Is it the day it is accepted for publication, the day an e-print is published online, or the day the physical paper copy is released? PubMed have thought about this and defined some rules -
http://www.nlm.nih.gov/bsd/licensee/elements_article_source.html
Essentially they are saying it is up to the publisher to decide when the date should be. An example is e.g.
http://www.ncbi.nlm.nih.gov/pubmed/20860432
The e-pub was in 2010, but the paper version was in 2011. The publishers (and hence PubMed) go for the later date in this case.
The pubdate value used in this pipeline is the most appropriate one based on what PubMed think it should be. This will likely be different to what is in the ALSPAC list of papers since stuff goes there as soon as it is available online.

It all depends on the \<Article PubModel="Print-Electronic"\> tag. This way round takes the print year, maybe we should use the e-pub year?

Some dates are even more complex and include a date range, see e.g. http://www.ncbi.nlm.nih.gov/pubmed/20648575 and http://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#articledate The only way I can see to automate processing of these is to pick the first four characters and hope it is a year!
