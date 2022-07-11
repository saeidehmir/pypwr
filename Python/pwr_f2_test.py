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

def pwr_f2_test(u=None, v=None, f2=None, sig_level=0.05, power=None):
    if sum(x is None for x in [u, v, f2, power, sig_level]) != 1:
        raise Exception("exactly one of u, v, f2, power, and sig.level must be NULL")
        sys.exit(0)
    if f2 is not None:
        if isinstance(f2,str):
            f2 = cohen_ES(test="f2",size=f2)["effect_size"]
        if f2<0:
            raise Exception("f2 must be positive")
            sys.exit(0)
    if (u is not None) and u<1:
        raise Exception("degree of freedom u for numerator must be at least 1")
        sys.exit(0)
    if (v is not None) and v<1:
        raise Exception("degree of freedom v for denominator must be at least 1")
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
    p_body = "ncf.sf(f_dist.isf(sig_level, dfn=u, dfd=v), dfn=u, dfd=v, nc=f2*(u+v+1))"
    if power is None:
        power = eval(p_body)
    elif u is None:
        loc = {}
        target_function_def = "def target_function(u):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        u = brentq(loc["target_function"], 1+1e-10, 100)
    elif v is None:
        loc = {}
        target_function_def = "def target_function(v):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        v = brentq(loc["target_function"], 1+1e-10, 1e+09)
    elif f2 is None:
        loc = {}
        target_function_def = "def target_function(f2):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        f2 = brentq(loc["target_function"], 1e-07, 1e+07)
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
    METHOD = "Multiple regression power calculation"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"u":u, "v":v, "f2":f2, "sig_level":sig_level, "power":power,
            "method":METHOD}
