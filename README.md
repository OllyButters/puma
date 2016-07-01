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


