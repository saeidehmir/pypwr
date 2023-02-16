# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import sys
import numpy as np
from cohen_ES import cohen_ES
from scipy.stats import norm
from scipy.optimize import brentq


def pwr_2p_test(alternative, h=None, n=None, sig_level=0.05, power=None):
    assert alternative in ["two.sided","less","greater"]
    alternative_options = ["less", "two.sided", "greater"]
    if sum(x is None for x in [h,n,power, sig_level]) != 1:
        raise Exception("exactly one of h, n, power, and sig_level must be NULL")
        sys.exit(0)
    if (h is not None) and isinstance(h,str):
        h = cohen_ES(test="p",size=h)["effect_size"]
    if (n is not None) and n<1:
        raise Exception("number of observations in each group must be at least 1")
        sys.exit(0)
    if (sig_level is not None):
        #print(type(sig_level))
        if type(sig_level)!=int and type(sig_level)!=float:
            raise Exception("sig_level should be numerical")
            sys.exit(0)
        elif sig_level<0 or sig_level>1:
            raise Exception("sig_level should be in [0, 1]")
            sys.exit(0)
    if (power is not None):
        if type(power)!=int and type(power)!=float:
            raise Exception("power should be numerical")
            sys.exit(0)
        elif power<0 or power>1:
            raise Exception("power should be in [0, 1]")
            sys.exit(0)
    tside = alternative_options.index(alternative)
    if tside==1 and (h is not None):
        h = np.abs(h)
    if tside == 1:
        p_body = "norm.sf(norm.isf(sig_level/2)-h*np.sqrt(n/2)) + norm.cdf(norm.ppf(sig_level/2)-h*np.sqrt(n/2))"
    elif tside == 2:
        p_body = "norm.sf(norm.isf(sig_level)-h*np.sqrt(n/2))"
    elif tside == 0:
        p_body = "norm.cdf(norm.ppf(sig_level)-h*np.sqrt(n/2))"
    if power is None:
        power = eval(p_body)
    elif h is None:
        loc = {}
        target_function_def = "def target_function(h):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        if tside==1:
            h = brentq(loc["target_function"], 1e-10, 10)
        elif tside==0:
            h = brentq(loc["target_function"], -10,5)
        elif tside==2:
            h = brentq(loc["target_function"], -5, 10)
    elif n is None:
        loc = {}
        target_function_def = "def target_function(n):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        n = brentq(loc["target_function"], 2+1e-10, 1e+09)
    elif sig_level is None:
        loc = {}
        target_function_def = "def target_function(sig_level):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        sig_level = brentq(loc["target_function"], 1e-10, 1-1e-10)
    else:
        raise Exception("internal error")
        sys.exit(0)
    NOTE = "same sample sizes"
    METHOD = "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"h":h, "n":n, "sig_level":sig_level, "power":power,
            "alternative":alternative, "method":METHOD, "note":NOTE}