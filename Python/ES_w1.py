# -*- coding: utf-8 -*-
"""
Created in Jul 2022

@author: Ernest Namdar
"""
import numpy as np

def ES_w1(P0,P1):
    return np.sqrt(np.sum((P1-P0)**2/P0))
