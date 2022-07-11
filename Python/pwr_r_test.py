# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import sys
import numpy as np
from cohen_ES import cohen_ES
from scipy.stats import t
from scipy.stats import norm
from scipy.optimize import brentq

def pwr_r_test(alternative, n=None, r=None, sig_level=0.05, power=None):
    assert alternative in ["two.sided","less","greater"]
    alternative_options = ["less", "two.sided", "greater"]
    if sum(x is None for x in [n, r, power, sig_level]) != 1:
        raise Exception("exactly one of n, r, power, and sig.level must be NULL")
        sys.exit(0)
    if (r is not None) and isinstance(r,str):
        r = cohen_ES(test="r",size=r)["effect_size"]
    if (n is not None) and n<4:
        raise Exception("number of observations must be at least 4")
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
    if tside==1 and (r is not None):
        r = np.abs(r)
    if tside == 1:
        p_body = "norm.cdf(((np.arctanh(r) + r/(2*(n-1)))-(np.arctanh(np.sqrt((t.isf(sig_level/2, df=n-2))**2/((t.isf(sig_level/2, df=n-2))**2+n-2)))))*np.sqrt(n-3))+norm.cdf((-(np.arctanh(r) + r/(2*(n-1)))-(np.arctanh(np.sqrt((t.isf(sig_level/2, df=n-2))**2/((t.isf(sig_level/2, df=n-2))**2+n-2)))))*np.sqrt(n-3))"
        # ttt = t.isf(sig_level/2, df=n-2)
        # rc = np.sqrt(ttt**2/(ttt**2+n-2))
        # zr = np.arctanh(r) + r/(2*(n-1))
        # zrc = np.arctanh(rc) # + rc/(2*(n-1))
        # p_body = norm.cdf((zr-zrc)*np.sqrt(n-3))+norm.cdf((-zr-zrc)*np.sqrt(n-3))

    elif tside == 2:
        p_body = "norm.cdf(((np.arctanh(r) + r/(2*(n-1)))-(np.arctanh(np.sqrt((t.isf(sig_level, df=n-2))**2/((t.isf(sig_level, df=n-2))**2+n-2)))))*np.sqrt(n-3))"
        # ttt = t.isf(sig_level, df=n-2)
        # rc = np.sqrt(ttt**2/(ttt**2+n-2))
        # zr = np.arctanh(r) + r/(2*(n-1))
        # zrc = np.arctanh(rc) # + rc/(2*(n-1))
        # p_body = norm.cdf((zr-zrc)*np.sqrt(n-3))

    elif tside == 0:
        p_body = "norm.cdf(((np.arctanh(-r) - r/(2*(n-1)))-(np.arctanh(np.sqrt((t.isf(sig_level, df=n-2))**2/((t.isf(sig_level, df=n-2))**2+n-2)))))*np.sqrt(n-3))"

        # r = -r
        # ttt = t.isf(sig_level, df=n-2)
        # rc = np.sqrt(ttt**2/(ttt**2+n-2))
        # zr = np.arctanh(r) + r/(2*(n-1))
        # zrc = np.arctanh(rc) # + rc/(2*(n-1))
        # p_body = norm.cdf((zr-zrc)*np.sqrt(n-3))

    if power is None:
        power = eval(p_body)
    elif r is None:
        loc = {}
        target_function_def = "def target_function(r):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        if tside==1:
            r = brentq(loc["target_function"], 1e-10, 1-1e-10)
        elif tside==0 or tside==2:
            r = brentq(loc["target_function"], -1+1e-10, 1-1e-10)
    elif n is None:
        loc = {}
        target_function_def = "def target_function(n):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        n = brentq(loc["target_function"], 4+1e-10, 1e+09)
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
    METHOD = "approximate correlation power calculation (arctangh transformation)"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"n":n, "r":r, "sig_level":sig_level, "power":power,
            "alternative":alternative, "method":METHOD}
