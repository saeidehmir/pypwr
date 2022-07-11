# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import numpy as np

def ES_w2(P):
    pi = np.sum(P, axis=0)
    pj = np.sum(P, axis=1)
    P0 = np.matmul(pi,pj)
    return np.sqrt(np.sum((P-P0)**2/P0))
