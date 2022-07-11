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

def pwr_norm_test(alternative, d=None, n=None, sig_level=0.05, power=None):
    assert alternative in ["two.sided","less","greater"]
    alternative_options = ["less", "two.sided", "greater"]
    if sum(x is None for x in [d, n, power, sig_level]) != 1:
        raise Exception("exactly one of d, n, power, and sig.level must be NULL")
        sys.exit(0)
    if (d is not None) and isinstance(d,str):
        d = cohen_ES(test="t",size=d)["effect_size"]
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
    if tside==1 and (d is not None):
        d = np.abs(d)
    if tside == 1:
        p_body = "norm.sf(norm.isf(sig_level/2)-d*np.sqrt(n)) + norm.cdf(norm.ppf(sig_level/2)-d*np.sqrt(n))"
    elif tside == 2:
        p_body = "norm.sf(norm.isf(sig_level)-d*np.sqrt(n))"
    elif tside == 0:
        p_body = "norm.cdf(norm.ppf(sig_level)-d*np.sqrt(n))"
    if power is None:
        power = eval(p_body)
    elif d is None:
        loc = {}
        target_function_def = "def target_function(d):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        if tside==1:
            d = brentq(loc["target_function"], 1e-10, 10)
        elif tside==0:
            d = brentq(loc["target_function"], -10,5)
        elif tside==2:
            d = brentq(loc["target_function"], -5, 10)
    elif n is None:
        loc = {}
        target_function_def = "def target_function(n):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        n = brentq(loc["target_function"], 1+1e-10, 1e+09)
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
    METHOD = "Mean power calculation for normal distribution with known variance"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"d":d, "n":n, "sig_level":sig_level, "power":power,
            "alternative":alternative, "method":METHOD}
