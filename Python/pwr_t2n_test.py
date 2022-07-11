# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import sys
import numpy as np
from cohen_ES import cohen_ES
from scipy.stats import t
from scipy.stats import nct
from scipy.stats import norm
from scipy.optimize import brentq

def pwr_t2n_test(alternative, n1=None, n2=None, d=None, sig_level=0.05, power=None):
    assert alternative in ["two.sided","less","greater"]
    alternative_options = ["less", "two.sided", "greater"]
    if sum(x is None for x in [n1, n2, d, power, sig_level]) != 1:
        raise Exception("exactly one of n1, n2, d, power, and sig.level must be NULL")
        sys.exit(0)
    if (d is not None) and isinstance(d,str):
        d = cohen_ES(test="t",size=d)["effect_size"]
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
    if (n1 is not None) and n1<2:
        raise Exception("number of observations in the first group must be at least 2")
        sys.exit(0)
    if (n2 is not None) and n2<2:
        raise Exception("number of observations in the second group must be at least 2")
        sys.exit(0)
    tsample = 2

    ttside = alternative_options.index(alternative)
    if ttside == 1:
        tside = 1
    else:
        tside=0
    tside += 1 #to avoid div by 0
    if tside==2 and (d is not None):
        d = np.abs(d)
    if ttside == 0:
        p_body = "nct.cdf(t.ppf(sig_level/tside, df=n1+n2-2), df=n1+n2-2, nc=d*(1/np.sqrt(1/n1+1/n2)))"
        # nu = n1+n2-2
        # qu = t.ppf(sig_level/tside, df=nu)
        # nct.cdf(qu, df=nu, nc=d*(1/np.sqrt(1/n1+1/n2)))
    elif ttside == 1:
        p_body = "nct.sf(t.isf(sig_level/tside, df=n1+n2-2), df=n1+n2-2, nc=d*(1/np.sqrt(1/n1+1/n2)))+nct.cdf(-(t.isf(sig_level/tside, df=n1+n2-2)), df=n1+n2-2, nc=d*(1/np.sqrt(1/n1+1/n2)))"
        # nu = n1+n2-2
        # qu = t.isf(sig_level/tside, df=nu)
        # nct.sf(qu, df=nu, nc=d*(1/np.sqrt(1/n1+1/n2))+nct.cdf(-qu, df=nu, nc=d*(1/np.sqrt(1/n1+1/n2))
    elif ttside == 2:
        p_body = "nct.sf(t.isf(sig_level/tside, df=n1+n2-2), df=n1+n2-2, nc=d*(1/np.sqrt(1/n1+1/n2)))"
        # nu = n1+n2-2
        # qu = t.isf(sig_level/tside, df=nu)
        # nct.sf(qu, df=nu, nc=d*(1/np.sqrt(1/n1+1/n2)))
    if power is None:
        power = eval(p_body)
    elif d is None:
        loc = {}
        target_function_def = "def target_function(d):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        if tside-1==0:
            d = brentq(loc["target_function"], -10, 5)
        elif tside-1==1:
            d = brentq(loc["target_function"], 1e-07, 10)
        elif tside-1==2:
            d = brentq(loc["target_function"], -5, 10)
    elif n1 is None:
        loc = {}
        target_function_def = "def target_function(n1):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        n1 = brentq(loc["target_function"], 2+1e-10, 1e+09)
    elif n2 is None:
        loc = {}
        target_function_def = "def target_function(n2):return "+p_body+"-power"
        exec(target_function_def, loc)
        #print(loc)
        loc.update(globals())
        loc.update(locals())
        n2 = brentq(loc["target_function"], 2+1e-10, 1e+09)
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
    METHOD = "t test power calculation"
    # print("h=",h, "n=", n, "power=",power, "sig_level=", sig_level)
    return {"n1":n1, "n2":n2, "d":d, "sig_level":sig_level, "power":power,
            "alternative":alternative, "method":METHOD}
