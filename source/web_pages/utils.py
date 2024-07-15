############################################################
# Sort an input list of papers by the sort parameter
############################################################
def sort_hashes_by(papers, hashes, sort_by):
    temp_papers = {}
    for this_paper in papers:
        if this_paper['IDs']['hash'] in hashes:
            if sort_by == 'year':
                temp_papers[this_paper['IDs']['hash']] = (this_paper['clean']['clean_date']['year'])
    sorted_hashes = [k for k, v in sorted(temp_papers.items(), key=lambda item: item[1], reverse=True)]
    return sorted_hashes

###########################################################
# Util to stick in commas between 1000s in big integers
###########################################################
def intWithCommas(x):
    if not isinstance(x, int):
        raise TypeError("Parameter must be an integer.")
    if x < 0:
        return '-' + intWithCommas(-x)
    result = ''
    while x >= 1000:
        x, r = divmod(x, 1000)
        result = ",%03d%s" % (r, result)
    return "%d%s" % (x, result)

