#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import csv

# Good plots
# NCDS
# keywords @ 0.7 has step change in it

# input_file = '../../data/alspac/title_lemmatized_by_year_weighted.csv'
# input_file = '../../data/ncds/title_lemmatized_by_year_weighted.csv'
# input_file = '../../data/alspac/keywords_lemmatized_by_year_weighted.csv'
input_file = '../../data/ncds/keywords_lemmatized_by_year_weighted.csv'
input_file = '../../data/ncds/temp.csv'
# input_file = '../../data/alspac/abstract_lemmatized_by_year_weighted.csv'
# input_file = '../../data/ncds/abstract_lemmatized_by_year_weighted.csv'

# I think was gets lemmatized to wa
ignore_list = ['and', 'of', 'in', 'a', 'the', 'study', 'with', 'at', 'to',
               'from', 'on', 'team', 'between', 'were', 'for', 'that', 'by',
               'n', 'y', 'human', 'study', 'be', '0', '1', 'or', 'six', 'will',
               '=', 'is', 'we', 'wa']


# minimum that the MAX value for a word must have.
# ncds
# title =  0.25/ keywords = 0.7 / abstrac = 1.5
freq_threshold = 0.7


# alspac
# title = 0.25 / keywords = 1 / abstrac = 1.5
# freq_threshold = 0.25


# Slice the data by these years. Min will be kept, max will be deleted
# ncds
min_year = '1996'
max_year = '2016'

# alspac
# min_year = '1996'
# max_year = '2016'

with open(input_file) as f:
    data = csv.reader(f)

    # The first row is the years
    years = next(data)
    years.pop(0)

    print(years)

    # Find the location  of the MIN year in the list
    min_year_location = years.index(min_year)
    print('Location of min year: ' + str(min_year_location))

    # Delete the data from the beginning of list to where the min year is
    years[:min_year_location] = []

    print(years)

    # Find the location  of the MAX year in the list
    max_year_location = years.index(max_year)
    print('Location of max year: ' + str(max_year_location))

    # Delete the data from the beginning of list to where the min year is
    years[max_year_location:] = []

    print(years)

    ok_freqs = 0
    # Cycle through each data line plotting the frequencies
    for row in data:
        freqs = row

        # Skip common words like 'the'
        word = row.pop(0)
        if word in ignore_list:
            continue

        # Delete the data from the beginning of list to where the min year is
        freqs[0:min_year_location] = []

        # Delete the data from the max years location to the end
        freqs[max_year_location:] = []

        # If no frequencies are above the threashold then skip it
        if float(max(freqs)) < freq_threshold:
            continue

        plt.plot(years, freqs, label=word)
        ok_freqs = ok_freqs + 1

print('Number of words plotted: ' + str(ok_freqs))

# Put a date range in. Good to ignore early and late years as low numbers of papers.
plt.legend(fontsize='x-large')
plt.xlim(int(min_year), int(max_year))
plt.xticks(np.arange(int(min_year), int(max_year)+1, 5.0), fontsize='x-large')
plt.xlabel('Year', fontsize='x-large')
plt.yticks(fontsize='x-large')
plt.ylabel('Weighted word freqency', fontsize='x-large')
plt.show()
