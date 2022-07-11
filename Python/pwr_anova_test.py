# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import sys
import numpy as np
from cohen_ES import cohen_ES
from scipy.stats import f as f_dist
from scipy.stats import ncf
from scipy.optimize import brentq

def pwr_anova_test(k=None, n=None, f=None, sig_level=0.05, power=None):
    if sum(x is None for x in [k, n, f, power, sig_level]) != 1:
        raise Exception("exactly one of k, n, f, power, and sig_level must be NULL")
        sys.exit(0)
    if f is not None:
        if isinstance(f,str):
            f = cohen_ES(test="anov",size=f)["effect_size"]
        if f<0:
            raise Exception("f must be positive")
            sys.exit(0)
    if (k is not None) and k<2:
        raise Exception("number of groups must be at least 2")
        sys.exit(0)
    if (n is not None) and n<2:
        raise Exception("number of observations in each group must be at least 2")
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
    p_body = "ncf.sf(f_dist.isf(sig_level, dfn=k-1, dfd=(n-1)*k), dfn=k-1, dfd=(n-1)*k, nc=k*n*(f**2))"
    if power is None:
        power = eval(p_body)
    elif k is None:
        loc = {}
        target_function_def = "def target_function(k):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        k = brentq(loc["target_function"], 2+1e-10, 100)
    elif n is None:
        loc = {}
        target_function_def = "def target_function(n):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        n = brentq(loc["target_function"], 2+1e-10, 1e+09)
    elif f is None:
        loc = {}
        target_function_def = "def target_function(f):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        f = brentq(loc["target_function"], 1e-07, 1e+07)
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
    NOTE = "n is number in each group"
    METHOD = "Balanced one-way analysis of variance power calculation"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"k":k, "n":n, "f":f, "sig_level":sig_level, "power":power,
            "method":METHOD, "note":NOTE}
