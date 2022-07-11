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
from pwr_t_test import pwr_t_test
from pwr_t2n_test import pwr_t2n_test
from pwr_2p_test import pwr_2p_test
from pwr_2p2n_test import pwr_2p2n_test
from pwr_anova_test import pwr_anova_test
from pwr_chisq_test import pwr_chisq_test
from pwr_norm_test import pwr_norm_test
from pwr_p_test import pwr_p_test
from pwr_r_test import pwr_r_test
import matplotlib.pyplot as plt

def plot_power_htest(x):
    # initial checks
    if type(x) != dict:
        raise Exception("argument must be of type dictionary")
        sys.exit(0)
    pwr_methods = ["One-sample t test power calculation",
                   "Two-sample t test power calculation",
                   "Paired t test power calculation",
                   "t test power calculation",
                   "Difference of proportion power calculation for binomial distribution (arcsine transformation)",
                   "difference of proportion power calculation for binomial distribution (arcsine transformation)",
                   "Balanced one-way analysis of variance power calculation",
                   "Chi squared power calculation",
                   "Mean power calculation for normal distribution with known variance",
                   "proportion power calculation for binomial distribution (arcsine transformation)",
                   "approximate correlation power calculation (arctangh transformation)"]
    if not x["method"] in pwr_methods:
        raise Exception("the method ", x["method"], " is not supported. Supported methods include:", pwr_methods)
        sys.exit(0)
    # settings
    breaks = 20


    # case: One-sample, Two-sample or Paired t test
    if x["method"] == "One-sample t test power calculation" or x["method"] == "Two-sample t test power calculation" or x["method"] == "Paired t test power calculation":
        if x["method"] == "One-sample t test power calculation":
            x["test_type"]="one.sample"
        if x["method"] == "Two-sample t test power calculation":
            x["test_type"]="two.sample"
        if x["method"] == "Paired t test power calculation":
            x["test_type"]="paired"
        n = x["n"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_t_test(x["alternative"], x["test_type"], n=sample_size, d=x["d"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = x["method"]
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\neffect size d =" + str(round(x["d"],4)) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(np.ceil(n)) + "\n" + str(x["note"])


    # case: Two-sample t test with n1 and n2
    elif x["method"] == "t test power calculation":
        n = x["n1"] + x["n2"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n
        n_rel = x["n1"] / n # relative sample size; will be kept constant in claculations

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            n1 = np.ceil(sample_size*n_rel)
            n2 = sample_size - n1
            if n1<2 or n2<2:
                data_power.append(None)
            else:
                data_power.append(pwr_t2n_test(alternative=x["alternative"], n1=n1, n2=n2, d=x["d"], sig_level=x["sig_level"], power=None)["power"])

        # create labels
        title_string = x["method"]
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\neffect size d =" + str(round(x["d"],4)) + "\nalpha =" + str(x["sig_level"]) + "\nn1/n2 = " + str(round(n_rel, 2))
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(x["n1"]) + '+' + str(x["n2"]) + '=' + str(n)


    # case: Difference of proportion (same sample size)
    elif x["method"] == "Difference of proportion power calculation for binomial distribution (arcsine transformation)":
        n = x["n"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_2p_test(x["alternative"], n=sample_size, h=x["h"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = "Difference of proportion power calculation\nfor binomial distribution (arcsine transformation)"
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\neffect size h =" + str(round(x["h"],4)) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(np.ceil(n)) + "\n" + str(x["note"])


    # case: difference of proportion (different sample size)
    elif x["method"] == "difference of proportion power calculation for binomial distribution (arcsine transformation)":
        n = x["n1"] + x["n2"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n
        n_rel = x["n1"] / n # relative sample size; will be kept constant in claculations

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            n1 = np.ceil(sample_size*n_rel)
            n2 = sample_size - n1
            if n1<2 or n2<2:
                data_power.append(None)
            else:
                data_power.append(pwr_2p2n_test(alternative=x["alternative"], n1=n1, n2=n2, h=x["h"], sig_level=x["sig_level"], power=None)["power"])

        # create labels
        title_string = "Difference of proportion power calculation\nfor binomial distribution (arcsine transformation)"
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\neffect size h =" + str(round(x["h"],4)) + "\nalpha =" + str(x["sig_level"]) + "\nn1/n2 = " + str(round(n_rel, 2))
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(x["n1"]) + '+' + str(x["n2"]) + '=' + str(n)


    # case: ANOVA
    elif x["method"] == "Balanced one-way analysis of variance power calculation":
        n = x["n"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_anova_test(n=sample_size, k=x["k"], f=x["f"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = "Balanced one-way analysis of variance \npower calculation"
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "groups k =" + str(x["k"]) + "\neffect size f =" + str(round(x["f"],4)) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(np.ceil(n)) + "\n" + str(x["note"])




    # case: Chi Squared
    elif x["method"] == "Chi squared power calculation":
        n = x["N"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_chisq_test(N=sample_size, w=x["w"], df=x["df"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = x["method"]
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "effect size w =" + str(round(x["w"],4)) + "\ndf =" + str(x["df"]) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nN = " + str(np.ceil(n)) + "\n" + str(x["note"])



    # case: Normal distribution
    elif x["method"] == "Mean power calculation for normal distribution with known variance":
        n = x["n"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_norm_test(alternative=x["alternative"], n=sample_size, d=x["d"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = "Mean power calculation for normal distribution\nwith known variance"
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\neffect size d =" + str(round(x["d"],4)) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(np.ceil(n))


    # case: proportion
    elif x["method"] == "proportion power calculation for binomial distribution (arcsine transformation)":
        n = x["n"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_p_test(alternative=x["alternative"], n=sample_size, h=x["h"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = "proportion power calculation\nfor binomial distribution (arcsine transformation)"
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\neffect size h =" + str(round(x["h"],4)) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(np.ceil(n))


    # case: correlation
    elif x["method"] == "approximate correlation power calculation (arctangh transformation)":
        n = x["n"]
        n_upper = max(n*1.5, n+30) # upper at least 30 above n

        # generate data
        sample_sizes = list(np.arange(10, n_upper, (n_upper-10)/breaks))
        data_sample_size = np.array(sample_sizes)
        data_power = []
        for sample_size in sample_sizes:
            data_power.append(pwr_r_test(alternative=x["alternative"], n=sample_size, r=x["r"], sig_level=x["sig_level"], power=None)["power"])
        data_power = np.array(data_power)
        #print(data_power)

        # create labels
        title_string = "approximate correlation power calculation\n(arctangh transformation)"
        #legend_string = '\n'.join(("tails =" + str(x["alternative"]), r"effect size d =", str(round(x["d"],4)), r"alpha =", str(x["sig_level"])))
        legend_string = "tails =" + str(x["alternative"]) + "\nr =" + str(round(x["r"],4)) + "\nalpha =" + str(x["sig_level"])
        xlab_string = "sample size"
        ylab_string = "test power = 1 - "+ r'$\beta$'
        optimal_string = "optimal sample size \nn = " + str(np.ceil(n))
    else:
        raise Exception("unexpected case")
        sys.exit(0)


    # plot
    fig, ax = plt.subplots()
    plt.plot(data_sample_size,data_power, c="red", marker='o', markerfacecolor="black", markeredgecolor="black")
    plt.axvline(x=np.ceil(n), c="blue", linestyle=':')
    plt.xlabel(xlab_string)
    plt.ylabel(ylab_string)
    plt.title(title_string)
    plt.text(0.1, 0.9, legend_string, fontsize=12, transform=ax.transAxes)
    plt.text(0.85, 0.1, optimal_string, fontsize=12, transform=ax.transAxes)
    plt.gca().set_yticklabels([f'{x:.0%}' for x in plt.gca().get_yticks()])
