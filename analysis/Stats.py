# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 22:13:43 2021

@author: Ian
"""


import matplotlib.pyplot as plt
import pdf2image
import math
import numpy as np
import seaborn as sea
import matplotlib as mpl
import statistics as stat
import scipy.stats as st
import os
import gc
from matplotlib.backends.backend_pgf import FigureCanvasPgf
mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pgf.texsystem'] = 'lualatex'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams.update({'figure.autolayout': True})
mpl.rcParams.update({'font.size': 10})


def pdfpng(path):  # convert all pdfs in directory to png
    files = os.listdir(path)
    filelist = [(path + "/" + i) for i in files]
    print("Converting PDF to PNG")
    for f in filelist:
        if "pdf" in f:
            pdfsplit = f.split(".pdf")[0]
            imgpath = pdfsplit.split("/")[-1]
            pdf2image.convert_from_path(f, dpi=500, output_file=imgpath,
                                        fmt="png", output_folder=path)


def meanTimePlot(k, time, y, tau, labels, a=1):  # Plot mean value against time
    # labels=[ylabel,property scale in latex form (e.g. $\tau_m=0.5 \tau_{d0}$)]
    # print("Plotting mean quantities...")
    yl = labels[0]
    prop = labels[1]
    axlabel = labels[2]
    titlestr = "Variation of mean droplet {} with time (normalised)".format(
        prop)
    # Normalise data
    ynorm = []
    tnorm = []
    y0 = y[0]  # initial data
    for i in range(0, len(time)):
        tnorm.append(time[i] / tau)
        ynorm.append(y[i] / y0)
    fig = plt.figure(k)
    plt.plot(tnorm, ynorm, label=axlabel, alpha=a)
    plt.xlabel(r"$\frac{t}{\tau_{d0}}$")
    plt.ylabel(yl)
    # plt.title(titlestr)
    plt.grid(True)
    plt.legend(prop={'size': 8})
    return fig


def sumPlot(k, time, y, tau, labels, a=1):  # Plot mean value against time
    # labels=[axlabel in latex form (e.g. $\tau_m=0.5 \tau_{d0}$)]
    # print("Plotting sum of droplet mass...")
    yl = r"$\sum \frac{m_d}{m_{d0}}$"
    axlabel = labels[0]
    titlestr = "Variation of sum droplet mass with time (normalised)"
    # Normalise data
    ynorm = []
    tnorm = []
    y0 = y[0]  # initial data
    for i in range(0, len(time)):
        tnorm.append(time[i] / tau)
        ynorm.append(y[i] / y0)
    fig = plt.figure(k)
    plt.plot(tnorm, ynorm, label=axlabel, alpha=a)
    plt.xlabel(r"$\frac{t}{\tau_{d0}}$")
    plt.ylabel(yl)
    # plt.title(titlestr)
    plt.grid(True)
    plt.legend()
    return fig


def meanScalePlot(scale, time, tau, labels):  # Plot mean value against scale
    # labels=[timeprop, category] (e.g. [evaporation time, mass timescale])
    timeprop = labels[0]
    category = labels[1]
    xl = labels[2]
    titlestr = "Variation of mean time taken to {} with {}".format(
        timeprop, category)
    tnorm = []  # normalise time data
    for i in time:
        tnorm.append(i / tau)
    plt.figure(dpi=500)
    plt.plot(scale, tnorm, "x-")
    plt.xlabel(xl)
    plt.ylabel(r"$\frac{t}{\tau{d0}}$")
    # plt.title(titlestr)
    plt.grid(True)
    plt.show()


def pdfPlot2d(dat, labels):  # 2D PDF, labels=[xlabel, title, filename]
    print("Plotting 2D PDF...")
    xl = labels[0]
    yl = "Density"
    tl = labels[1]
    filename = labels[2]
    var = stat.variance(dat)
    sk = st.skew(dat)
    pdf2d = plt.figure(dpi=250)
    sea.displot(dat, kind="kde")
    plt.xlabel(xl)
    # plt.title(tl)
    plt.grid(True)
    plt.annotate("Variance = {:.4e}\nSkewness = {:.4e}".format(
        var, sk), (0.05, 0.9), xycoords="axes fraction")
    try:
        plt.savefig(filename)
        plt.close(pdf2d)
    except ValueError:
        print("Savefig failed:", filename)
        pass


def pdfPlot3d(xin, yin, labels):  # 3D PDF, labels=[xlabel, ylabel, title, filename]
    # try:
    print("Plotting 3D PDF...")
    xl = labels[0]
    yl = labels[2]
    tl = labels[1]
    filename = labels[3]
    xvar = stat.variance(xin)
    xsk = st.skew(xin)
    yvar = stat.variance(yin)
    ysk = st.skew(yin)
    pdf3d = plt.figure(dpi=250)
    # sea.kdeplot(x=xin, y=yin, bw_adjust=1.75)
    sea.kdeplot(x=xin, y=yin)
    plt.grid(True)
    plt.xlabel(xl)
    plt.ylabel(yl)
    # plt.title(tl)
    plt.annotate("Variance = ({:.4e}, {:.4e})\nSkewness = ({:.4e}, {:.4e})".format(xvar, yvar, xsk, ysk),
                 (0.05, 0.9), xycoords="axes fraction")
    try:
        plt.savefig(filename)
        plt.close(pdf3d)
    except ValueError:
        print("Savefig failed:", filename)
        # pdf3d.show()
        pass


def newdir(directory):
    try:
        os.mkdir(directory)
    except OSError:
        pass


# obtain all data directories
simdir_list = [name+"/" for name in os.listdir()]
simdir_list.remove("stokes_numbers/")  # irrelevant sim to this IP
simdir_list.remove("m1m2_perf/")  # obtain data separately
simdir_list.remove("Stats.py/")
dirlist = []
for i in simdir_list:
    simtypes = os.listdir(i)
    if "data" in simtypes:
        dirstr = i + "data/"
        dirlist.append(dirstr)
    else:
        for j in simtypes:
            dirstr = i + j + "/data/"
            dirlist.append(dirstr)
# dirlist.remove("timescales/tau_m/data/")
dirlist.remove("timescales/tau_v/data/")
dirlist.remove("timescales/tau_h/data/")
dirlist.remove("fluidprops/mu_G/data/")
dirlist.remove("fluidprops/rho_G/data/")
dirlist.remove("fluidprops/T_G/data/")
dirlist.remove("tgvmag/data/")

# obtain full data directory for each simulation scale
fulldir = []
for i in dirlist:
    scalelist = os.listdir(i)
    for j in scalelist:
        dirstr = i + j + "/"
        fulldir.append(dirstr)


# obtain text files
reset = False
scalelist = []
mevaplist = []
mtsslist = []
mv0list = []
alldata = {}
count = 1
alldir = len(fulldir)
for i in fulldir:
    directory = i
    # obtain string for plot title
    if "tgvmag" in i:
        chmtransfer = i.split("/")[2]
    else:
        chmtransfer = i.split("/")[3]
    uscore = chmtransfer.split("_")[-2:]
    scale = eval(".".join(uscore))
    scalelist.append(scale)
    if scale == 5.0:  # reset scalelist, mean lists at the end of loop
        reset = True
    if "tau_m" in i:
        propstr = r"$\tau_m = {} ".format(scale) + r"\tau_{d}$"
    elif "tau_h" in i:
        propstr = r"$\tau_h = {} ".format(scale) + r"\tau_{d}$"
    elif "tau_v" in i:
        propstr = r"$\tau_v = {} ".format(scale) + r"\tau_{d}$"
    elif "mu_G" in i:
        propstr = r"$\mu_G = {} ".format(scale) + r"\mu_{G0}$"
    elif "rho_G" in i:
        propstr = r"$\rho_G = {} ".format(scale) + r"\rho_{G0}$"
    elif "T_G" in i:
        propstr = r"$T_G = {} ".format(scale) + r"T_{G0}$"
    elif "tgv" in i:
        propstr = r"$M = {}$".format(scale)

    # prepare ylabel strings
    diameter = r"$\frac{D^2}{{D_0}^2}$"
    temperature = r"$\frac{{T_d}}{{T_{d0}}}$"
    velocity = r"$\frac{|v|}{|v_0|}$"

    # sort files and initialilse lists
    filelist = os.listdir(i)
    filedir = [i + file for file in filelist]

    setupfile = filedir[-1]  # obtain setup from setup file
    f = open(setupfile, "r")
    flines = f.readlines()
    f.close()
    taustr = flines[-1]
    taustr1 = taustr.split(" ")[-1]
    tau = float(taustr1.split("\n")[0])

    filedir.remove(filedir[-1])  # remove setup file after obtaining setup
    filedir.sort(key=os.path.getctime)
    tlist = []
    vabslist = []
    templist = []
    dsqlist = []
    masslist = []
    fvlist = []
    ftlist = []
    fdlist = []
    fmlist = []
    meanv = []
    meant = []
    meand = []
    meanm = []
    summ = []

    # Initialise plots
    dplot = plt.figure(0, dpi=500)
    tplot = plt.figure(1, dpi=500)
    vplot = plt.figure(2, dpi=500)
    splot = plt.figure(3, dpi=500)

    currfilenum = 1
    print("{} -- Reading files...".format(directory))
    for f in filedir:  # loop through each timestep file
        a = open(f, "r")
        lines = a.readlines()
        timestr = lines[0]
        time = float(timestr.split(" ")[-2])
        tlist.append(time)
        # print("{} -- T = {} -- Reading files... ({} %)".format(directory, time,
        #       round(currfilenum/len(filedir)*100, 3)), flush=True)
        # print("Time =", time)
        particles_list = [j for j in lines[1:-1]]
        a.close()
        vabs = []
        temp = []
        dsquare = []
        mass = []
        fullvabs = []
        fulltemp = []
        fulldsquare = []
        fullmass = []

        for particle in particles_list:
            # loop through each particle line in single file, split strings into numbers, ...
            # then put data into corresponding categories
            splitnewline = particle.split("\n")[0]
            splitcomma = splitnewline.split(",")
            vx = float(splitcomma[0])
            vy = float(splitcomma[1])
            vz = float(splitcomma[2])
            fvx = float(splitcomma[3])
            fvy = float(splitcomma[4])
            fvz = float(splitcomma[5])
            fullvabs.append(np.sqrt(((vx-fvx)**2)+((vy-fvy)**2)+((vz-fvz)**2)))
            fulltemp.append(float(splitcomma[6]))
            fulldsquare.append(float(splitcomma[7]))
            fullmass.append(float(splitcomma[8]))

            # filter out zero mass particles for mean calculation
            if float(splitcomma[8]) != 0:
                # velx.append(vx)
                # vely.append(vy)
                # velz.append(vz)
                vabs.append(np.sqrt(((vx-fvx)**2)+((vy-fvy)**2)+((vz-fvz)**2)))
                # fvelx.append(float(splitcomma[3]))
                # fvely.append(float(splitcomma[4]))
                # fvelz.append(float(splitcomma[5]))
                temp.append(float(splitcomma[6]))
                dsquare.append(float(splitcomma[7]))
                mass.append(float(splitcomma[8]))
        # mean properties
        try:
            dsqmean = stat.mean(dsquare)
        except stat.StatisticsError:
            dsqmean = meand[-1]
        try:
            tempmean = stat.mean(temp)
        except stat.StatisticsError:
            tempmean = meant[-1]
        try:
            massmean = stat.mean(mass)
        except stat.StatisticsError:
            massmean = 0
        try:
            vabsmean = stat.mean(vabs)
        except stat.StatisticsError:
            vabsmean = meanv[-1]

        # mass sum
        masssum = sum(mass)

        # append to full list (all timesteps)
        vabslist.append(vabs)
        templist.append(temp)
        dsqlist.append(dsquare)
        masslist.append(mass)

        fvlist.append(fullvabs)
        ftlist.append(fulltemp)
        fdlist.append(fulldsquare)
        fmlist.append(fullmass)

        meanv.append(vabsmean)
        meant.append(tempmean)
        meand.append(dsqmean)
        meanm.append(massmean)
        summ.append(masssum)

        currfilenum += 1

    print("{} -- Sorting particle data...".format(directory), flush=True)

    # SORT DATA INTO CHRONOLOGICAL ORDER
    tlist, vabslist, templist, dsqlist, masslist, fvlist, ftlist, fdlist, fmlist, meanv, meant, meand, meanm = zip(*sorted(zip(tlist, vabslist, templist, dsqlist, masslist, fvlist, ftlist,
                                                                                                                               fdlist, fmlist, meanv, meant, meand, meanm)))
    # dictionary of evap time + tss time + zero vel time of each particle, useful for scatter plots
    part = {}  # format {id1: [evap, tss, vel], id2:...}

    # determine evaporation time of each particle
    evaptime = []
    numpart = len(fmlist[0])
    # convert to list for easier operation of item removal
    partid = list(np.arange(0, numpart, 1))
    for n in partid:
        part[n] = [None, None, None]  # initialise dictionary

    for i in range(0, len(fmlist)):
        tpart = tlist[i]
        pmass = fmlist[i]  # list of particle mass of one timestep
        for j in partid:
            if pmass[j] == 0:  # if zero mass is achieved
                # remove particle id from array for next timestep
                partid.remove(j)
                evaptime.append(tpart)
                part[j][0] = tpart
    meanevap = stat.mean(evaptime)
    mevaplist.append(meanevap)

    # determine steady state temperature itme of each particle
    tsstime = []
    partid = list(np.arange(0, numpart, 1))  # must include to reset the array
    tss = meant[-1]  # steady state temperature
    for i in range(0, len(ftlist)):
        tpart = tlist[i]
        ptemp = ftlist[i]  # list of particle temp of one timestep
        pm = fmlist[i]
        for j in partid:
            ratio = abs(ptemp[j] - tss) / tss
            if ratio < 1e-2 and pm != 0:  # if tss is achieved
                # remove particle id from array for next timestep
                partid.remove(j)
                tsstime.append(tpart)
                part[j][1] = tpart
    meantss = stat.mean(tsstime)
    mtsslist.append(meantss)

    # determine zero relative velocity time of each particle
    v0time = []
    partid = list(np.arange(0, numpart, 1))  # must include to reset the array
    for i in range(0, len(fvlist)):
        tpart = tlist[i]
        pvel = fvlist[i]  # list of particle |v| of one timestep
        pm = fmlist[i]
        for j in partid:
            if pvel[j] == 0 and pm[j] != 0:
                partid.remove(j)
                v0time.append(tpart)
                part[j][2] = tpart
    if len(v0time) != 0:
        meanv0 = stat.mean(v0time)
    else:
        meanv0 = None
    mv0list.append(meanv0)

    print("{} -- Plotting...".format(directory))
    # Conditionals below are used to define types of plot to produce
    alpha = 1
    if "timescales" or "fluidprops" or "tgvmag" in f:
        if "tau_m" in f:
            dirnew = "timescales/tau_m/figures/"
            newdir(dirnew)
            pdfstr = "timescales/tau_m/figures/" + "taum_" + str(scale) + "_"
            scplotstr = "Mass timescale"
            scplotxl = r"Timescale coefficient $C_{\tau}$"
        elif "tau_h" in f:
            dirnew = "timescales/tau_h/figures/"
            newdir(dirnew)
            pdfstr = "timescales/tau_h/figures/" + "tauh_" + str(scale) + "_"
            scplotstr = "Temperature timescale"
            scplotxl = r"Timescale coefficient $C_{\tau}$"
            alpha = 0.5
        elif "tau_v" in f:
            dirnew = "timescales/tau_v/figures/"
            newdir(dirnew)
            pdfstr = "timescales/tau_v/figures/" + "tauv_" + str(scale) + "_"
            scplotstr = "Velocity timescale"
            scplotxl = r"Timescale coefficient $C_{\tau}$"
            # alpha = 0.5
        elif "mu_G" in f:
            dirnew = "fluidprops/mu_G/figures/"
            newdir(dirnew)
            pdfstr = "fluidprops/mu_G/figures/" + "mu_" + str(scale) + "_"
            scplotstr = "Fluid viscosity"
            scplotxl = r"Viscosity coefficient $C_{\mu}$"
            tau = scale * tau  # fix tau to reference viscosity
        elif "rho_G" in f:
            dirnew = "fluidprops/rho_G/figures/"
            newdir(dirnew)
            pdfstr = "fluidprops/rho_G/figures/" + "rho_" + str(scale) + "_"
            scplotstr = "Fluid density"
            scplotxl = r"Density coefficient $C_{\rho}$"
        elif "T_G" in f:
            dirnew = "fluidprops/T_G/figures/"
            newdir(dirnew)
            pdfstr = "fluidprops/T_G/figures/" + "tg_" + str(scale) + "_"
            scplotstr = "Fluid temperature"
            scplotxl = r"Temperature coefficient $C_{T}$"
            mu_G = 6.109e-6 + 4.604e-8 * 1000 - 1.051e-11 * \
                (1000**2)  # force reference viscosity, T_G = 1000K
            tau = scale * tau  # fix tau to reference temperature
        elif "tgvmag" in f:
            dirnew = "tgvmag/figures/"
            newdir(dirnew)
            pdfstr = "tgvmag/figures/" + "tgv_" + str(scale) + "_"
            scplotstr = "Taylor-Green Vortex magnitude"
            scplotxl = r"Magnitude $M$"

        # rearrange lists in chronological order
        tlist, meand, meant, meanv = zip(*sorted(zip(tlist, meand,
                                                     meant, meanv)))

        # mean droplet property plots (this may fail if code is also plotting PDFs)
        dplot = meanTimePlot(0, tlist, meand, tau, [
            diameter, "diameter", propstr], alpha)
        tplot = meanTimePlot(1, tlist, meant, tau, [
            temperature, "temperature", propstr], alpha)
        vplot = meanTimePlot(2, tlist, meanv, tau, [
            velocity, "absolute relative velocity", propstr],
            alpha)
        splot = sumPlot(3, tlist, summ, tau, [propstr], alpha)

        # mean time vs scale plots (this may fail if code is also plotting PDFs)
        if scale == 5.0:
            print(scalelist, mevaplist)
            alldata[scplotstr] = (scalelist, mevaplist)
            meanScalePlot(scalelist, mevaplist, tau, [
                          "evaporate droplets", scplotstr, scplotxl])
            meanScalePlot(scalelist, mtsslist, tau, [
                          "attain steady state temperature", scplotstr, scplotxl])

        # 3D PDF for some timesteps
        numtfiles = len(tlist)  # find number of files in this timescale
        # maxt = tlist(len(tlist)-1)
        numgraphs = 5  # number of PDFs to generate per timescale
        step = math.floor(numtfiles / numgraphs)
        # print(numtfiles, step)

        if scale == 0.5 or scale == 3 or scale == 5:
            pass
            # for i in range(0, numgraphs+1):
            #     if i == 0:
            #         ind = 1
            #     else:
            #         ind = step * i
            #     if i == 0:
            #         titlestr = r"$t \approx 0$"
            #     elif i == 5:
            #         titlestr = r"$t \approx t_{evap}$"
            #     else:
            #         titlestr = r"$t = {}\,$".format(
            #             str(round(1/numgraphs * (i), 2))) + r"$t_{evap}$"
            #     while ind >= len(tlist):
            #         ind -= 1
            #         titlestr = r"$t \approx t_{evap}$"
            #     while len(templist[ind]) <= 2:
            #         ind -= 1
            #         titlestr = r"$t \approx t_{evap}$"
            #     temp = templist[ind]
            #     dsq = dsqlist[ind]
            #     tnorm = [(j / meant[-1]) for j in temp]
            #     dnorm = [(j / meand[0]) for j in dsq]
            #     filename = pdfstr + "d2td_" + str(i) + ".pdf"
            #     print(i, ind, numtfiles, filename)
            #     pdfPlot3d(tnorm, dnorm, [r"$\frac{T}{T_{ss}}$", titlestr,
            #                              r"$\frac{D^2}{{D_0}^2}$", filename])

            # for i in range(0, numgraphs+1):
            #     if i == 0:
            #         ind = 1
            #     else:
            #         ind = step * i
            #     if i == 0:
            #         titlestr = r"$t \approx 0$"
            #     elif i == 5:
            #         titlestr = r"$t \approx t_{evap}$"
            #     else:
            #         titlestr = r"$t = {}\,$".format(
            #             str(round(1/numgraphs * (i), 2))) + r"$t_{evap}$"
            #     while ind >= len(tlist):
            #         ind -= 1
            #         titlestr = r"$t \approx t_{evap}$"
            #     while len(vabslist[ind]) <= 2:
            #         ind -= 1
            #         titlestr = r"$t \approx t_{evap}$"
            #     vabs = vabslist[ind]
            #     dsq = dsqlist[ind]
            #     v0 = vabslist[0][0]
            #     vnorm = [(j/v0) for j in vabs]
            #     dnorm = [(j / meand[0]) for j in dsq]
            #     filename = pdfstr + "vabs_" + str(i) + ".pdf"
            #     print(i, ind, numtfiles, filename)
            #     pdfPlot3d(vnorm, dnorm, [r"$\frac{|v|}{|v_0|}$", titlestr,
            #                              r"$\frac{D^2}{{D_0}^2}$", filename])

        # 3D PDF on evap time-tss time correlation
            # pevt = []
            # ptsst = []
            # for i in part:
            #     p = part[i]
            #     evap = p[0]
            #     tss = p[1]
            #     try:
            #         pevnorm = evap / tau
            #         ptssnorm = tss / tau
            #         pevt.append(pevnorm)
            #         ptsst.append(ptssnorm)
            #     except TypeError:
            #         pass
            # filename = pdfstr + "evaptss" + ".pdf"
            # print(filename)
            # pdfPlot3d(ptsst, pevt, [r"$\frac{t_{ss}}{\tau}$", None,
            #                         r"$\frac{t_{evap}}{\tau}$", filename])

    if reset:
        scalelist = []
        mevaplist = []
        mtsslist = []
        mv0list = []
        reset = False
    count += 1

# all plots
plt.figure(dpi=1000)
for i in alldata:
    tnorm = []  # normalise time data
    s = alldata.get(i)[0]
    ty = alldata.get(i)[1]
    for j in ty:
        tnorm.append(j / tau)
    plt.plot(s, tnorm, "x-", label=i, linewidth=1)
plt.xlabel("Scale")
plt.ylabel(r"$\frac{t}{\tau{d0}}$")
plt.legend(prop={'size': 8})
plt.grid(True)
plt.show()


# m1m2_perf
dirname = "m1m2_perf/data/"
tfile = dirname + "times2.txt"
f = open(tfile, "r")
flines = f.readlines()
f.close()
n1list = []
t1list = []
ts1list = []
ls1list = []
n2list = []
t2list = []
ts2list = []
ls2list = []

for line in flines[1:]:
    a = line.split("\n")[0]
    model = int(a.split(",")[0])
    tstep = float(a.split(",")[1])
    lstep = int(a.split(",")[2])
    numpart = int(a.split(",")[3])
    t = float(a.split(",")[-1])
    # print(model, tstep, lstep, numpart, t)
    if model == 1:
        n1list.append(numpart)
        t1list.append(t)
        ts1list.append(tstep)
        ls1list.append(lstep)
    elif model == 2:
        n2list.append(numpart)
        t2list.append(t)
        ts2list.append(tstep)
        ls2list.append(lstep)

plt.figure(dpi=500)
plt.loglog(n1list, t1list, "x-", label="M1")
plt.loglog(n2list, t2list, "x-", label="M2")
plt.xlabel("Number of droplets")
plt.ylabel("Time")
plt.legend()
plt.grid(True, which="major")
plt.grid(True, which="minor", color="0.9")
plt.savefig("m1m2_perf/data/SimTime.png")
plt.show()

w = 1
plt.figure(dpi=500)
plt.loglog(n1list, t1list, "x--", label=r"M1, Logstep $=100\Delta t$", linewidth=w)
plt.loglog(n2list, t2list, "x-", label=r"M2, Logstep $=100\Delta t$", linewidth=w)
tfile = dirname + "times500.txt"
f = open(tfile, "r")
flines = f.readlines()
f.close()
n1list = []
t1list = []
ts1list = []
ls1list = []
n2list = []
t2list = []
ts2list = []
ls2list = []
for line in flines[1:]:
    a = line.split("\n")[0]
    model = int(a.split(",")[0])
    tstep = float(a.split(",")[1])
    lstep = int(a.split(",")[2])
    numpart = int(a.split(",")[3])
    t = float(a.split(",")[-1])
    if model == 1:
        n1list.append(numpart)
        t1list.append(t)
        ts1list.append(tstep)
        ls1list.append(lstep)
    elif model == 2:
        n2list.append(numpart)
        t2list.append(t)
        ts2list.append(tstep)
        ls2list.append(lstep)
plt.loglog(n1list, t1list, "x--", label=r"M1, Logstep $=500\Delta t$", linewidth=w)
plt.loglog(n2list, t2list, "x-", label=r"M2, Logstep $=500\Delta t$", linewidth=w)
tfile = dirname + "times1000.txt"
f = open(tfile, "r")
flines = f.readlines()
f.close()
n1list = []
t1list = []
ts1list = []
ls1list = []
n2list = []
t2list = []
ts2list = []
ls2list = []
for line in flines[1:]:
    a = line.split("\n")[0]
    model = int(a.split(",")[0])
    tstep = float(a.split(",")[1])
    lstep = int(a.split(",")[2])
    numpart = int(a.split(",")[3])
    t = float(a.split(",")[-1])
    if model == 1:
        n1list.append(numpart)
        t1list.append(t)
        ts1list.append(tstep)
        ls1list.append(lstep)
    elif model == 2:
        n2list.append(numpart)
        t2list.append(t)
        ts2list.append(tstep)
        ls2list.append(lstep)
plt.loglog(n1list, t1list, "x--", label=r"M1, Logstep $=1000\Delta t$", linewidth=w)
plt.loglog(n2list, t2list, "x-", label=r"M2, Logstep $=1000\Delta t$", linewidth=w)
plt.xlabel("Number of droplets")
plt.ylabel("Time")
plt.legend()
plt.grid(True, which="major")
plt.grid(True, which="minor", color="0.9")
plt.savefig("m1m2_perf/data/LogTime.png")
plt.show()

gc.collect()  # clear memory
print("\n**Reminder: Remember to remove plots regularly to save memory**")
