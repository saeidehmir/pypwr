# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
def cohen_ES(test, size):
    """

    Parameters
    ----------
    test : str
        Available options are "p","t","r","anov","chisq","f2".
    size : str
        Available options are "small","medium","large".

    Returns
    -------
    dict
        a dictionary of test, size, effect_size (calculated), and method.

    """
    assert test in ["p","t","r","anov","chisq","f2"]
    assert size in ["small","medium","large"]
    test_options = ["p","t","r","anov","chisq","f2"]
    # size_options = ["small","medium","large"]
    ntest = test_options.index(test)
    if ntest == 0:
        size_mapping = {"small":0.2, "medium":0.5, "large":0.8}
        ES = size_mapping[size]
    elif ntest == 1:
        size_mapping = {"small":0.2, "medium":0.5, "large":0.8}
        ES = size_mapping[size]
    elif ntest == 2:
        size_mapping = {"small":0.1, "medium":0.3, "large":0.5}
        ES = size_mapping[size]
    elif ntest == 3:
        size_mapping = {"small":0.1, "medium":0.25, "large":0.4}
        ES = size_mapping[size]
    elif ntest == 4:
        size_mapping = {"small":0.1, "medium":0.3, "large":0.5}
        ES = size_mapping[size]
    elif ntest == 5:
        size_mapping = {"small":0.02, "medium":0.15, "large":0.35}
        ES = size_mapping[size]
    METHOD = "Conventional effect size from Cohen (1982)"
    return {"test":test, "size":size, "effect_size":ES, "method":METHOD}

