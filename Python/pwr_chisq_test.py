# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import sys
import numpy as np
from cohen_ES import cohen_ES
from scipy.stats import chi2
from scipy.stats import ncx2
from scipy.optimize import brentq

def pwr_chisq_test(w=None, N=None, df=None, sig_level=0.05, power=None):
    if sum(x is None for x in [w, N, df, power, sig_level]) != 1:
        raise Exception("exactly one of w, N, df, power, and sig_level must be NULL")
        sys.exit(0)
    if w is not None:
        if isinstance(w,str):
            w = cohen_ES(test="chisq",size=w)["effect_size"]
        if w<0:
            raise Exception("w must be positive")
            sys.exit(0)
    if (N is not None) and N<1:
        raise Exception("number of observations must be at least 1")
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
    p_body = "ncx2.sf(chi2.isf(sig_level, df=df), df=df, nc=N*(w**2))"
    if power is None:
        power = eval(p_body)
    elif w is None:
        loc = {}
        target_function_def = "def target_function(w):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        w = brentq(loc["target_function"], 1e-10, 1e+09)
    elif N is None:
        loc = {}
        target_function_def = "def target_function(N):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        N = brentq(loc["target_function"], 1+1e-10, 1e+09)
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
    NOTE = "N is the number of observations"
    METHOD = "Chi squared power calculation"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"w":w, "N":N, "df":df, "sig_level":sig_level, "power":power,
            "method":METHOD, "note":NOTE}
