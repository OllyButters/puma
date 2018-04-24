#! /usr/bin/env python3

import matplotlib.pyplot as plt
import csv

# Could put a minumum threashold in this so if all values in list are below X then
# ignore it. Useful for when e.g. only a couple of mentions of word ever.

# input_file = '../data/alspac/title_lemmatized_by_year.csv'
# input_file = '../data/alspac/title_lemmatized_by_year_weighted.csv'
# input_file = '../data/alspac/keywords_lemmatized_by_year.csv'
# input_file = '../data/alspac/keywords_lemmatized_by_year_weighted.csv'
# input_file = '../data/alspac/abstract_lemmatized_by_year.csv'
input_file = '../data/ncds/abstract_lemmatized_by_year_weighted.csv'

ignore_list = ['and', 'of', 'in', 'a', 'the', 'study', 'with', 'at', 'to', 'from', 'on', 'team', 'between', 'were', 'for', 'that', 'by', 'n', 'y', 'human', 'study']

freq_threshold = 0.1

with open(input_file) as f:
    data = csv.reader(f)

    # The first row is the years
    years = next(data)
    years.pop(0)

    ok_freqs = 0
    # Cycle through each data line plotting the frequencies
    for row in data:
        freqs = row

        # Skip common words like the
        word = row.pop(0)
        if word in ignore_list:
            continue

        # If no frequencies are above the threashold then skip it
        if float(max(freqs)) < freq_threshold:
            continue

        plt.plot(years, freqs, label='bob')
        ok_freqs = ok_freqs + 1

print(ok_freqs)

# Put a date range in. Good to ignore early and late years as low numbers of papers.
# plt.xlim(1997, 2015)
plt.xlim(1970, 2015)

# Build the plot
plt.show()
