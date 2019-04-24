def count_ones(word_ones):
    (word, ones) = word_ones
    return word, sum(ones)


def process_entries(sample_list):
    """ Values that come out is not a list of dicts, as expected, but something called a Unwindowed value.
    Convert to list first, then looks like a list of dicts. So basically, you get a (key, [{dict1}, {dict2}, {dict3}])
    """
    (word, values) = sample_list
    list_dicts = list(values)

    return word, list_dicts
