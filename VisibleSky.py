# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 12:43:19 2020

@author: themi
"""

#######################################################################  
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
########################### WHAT IS THIS? #############################
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
#######################################################################
#                                                                     #
# This is a widget based off a script known as Revolving Door. So     #
# what does it do, you might ask? Select an observatory, any of CTIO, #
# VRO, or KNPO, and input a date to see the RA range of the telescope #
# or an set of RA and DEC coordinates to see if that point is visible #
# currently.                                                          #
#                                                                     #
#######################################################################
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
#######################################################################

import datetime
from astropy.time import Time
import numpy as np
import calendar
from calendar import monthrange
import csv
from string import whitespace
import pandas as pd
import argparse
import os.path
import sqlite3
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import math
import astroplan
import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TKAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from astroplan import download_IERS_A
download_IERS_A()

#######################################################################  
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
############################# FUNCTIONS ###############################
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
#######################################################################

#######################################################################
"""
*** What is this? ***
This is a set of functions which calculates the Greenwich Mean 
Time (GMST) given the year, month, day, hour, minute, and second.

*** How does it work? ***
The core function of this program is 'getGMST,' which takes 
utilizes several secondary functions to calculate the given GMST.
The secondary functions are J0, which calculates 'Julian day number 
at 0 h UT;' T0, which calculates the time in 'in Julian centuries 
between the Julian day J0 and J2000;' thetaG0, which 'Greenwich 
sidereal time [...] at 0 h UT;' and thetaG, which calculates 
'Greenwich sidereal time [...] at any other UT.'

*** Where are those quotations from and these formulae? ***
Right here: https://www.sciencedirect.com/book/9780081021330/orbital
-mechanics-for-engineering-students
"""

def J0(y,m,d): # where y is the year, m is the month, and d is the day
    
    int1 = int((m+9)/12)
    int2 = int(275*m/9)
    
    J0 = 367*y - int((7*y+int1)/4) + int2 + d + 1721013.5
    
    return J0

def T0(j0):
    
    T0 = (j0 - 2451545)/36525
    
    return T0

def thetaG0(t0):
    
    thetaG0 = 100.4606184 + 36000.77004*t0 + 0.000387933*t0**2 - ((2.58310)**(-8))*t0**3
    
    if thetaG0 < 0 or thetaG0 > 360:
        
        multiple = int(thetaG0/360)
        thetaG0 = thetaG0 - multiple*360
    
    return thetaG0

def thetaG(thetag0,h,minute,sec):
    
    thetaG = thetag0 + 360.98564724*(h + minute/60 + sec/3600)/24
    
    return thetaG

# For thetaG in hours, divide by 15
    
def getGMST(y,m,d,h,minute,sec): # For the gmst at the start of each month,
    #set day to 1, minute to 0, and second to 0
    
    y = float(y)
    m = float(m)
    d = float(d)
    minute = float(minute)
    sec = float(sec)
    
    thetag = thetaG(thetaG0(T0(J0(y,m,d))),h,minute,sec)
    GMST = thetag/15
    
    return GMST
#######################################################################

def getMST(gmst, Day, Hour, Min, Sec): # Determines means sidereal time in UTC
        
    GMST = gmst + (0.06571*float(Day)) + (1.002738*(float(Hour) + (float(Min)/60.0) + (float(Sec)/3600)))
    mst = (((GMST*15) - 70.80639)/15)%24
    
    return mst

def getDateTimefromJulian(julianDate):
    # Converts Julian date to UTC
                                                                                 
    con = sqlite3.connect(":memory:")
    query = "select datetime('%s')" %julianDate
    date = str(list(con.execute(query))[0][0])
   
    return date # UTC

def getNightRange(mstSet, mstRise): # East is the least for RA ???? 
    # Determines RA and DECs visible during given night
    
    if mstRise < mstSet:
        vRWest = mstSet
        vREast = mstRise + 24
    else:
        vRWest = mstSet
        vREast = mstRise

    visRange = (vRWest-5.25,vREast+5.25) # Max westbound vis, max eastbound viso
    
    v1 = []
    v2 = []
    v3 = []
    
    if visRange[0] <= 0 and visRange[1] >= 24:
        v2 = [0,24]
    elif visRange[0] < 0 and visRange[1] < 24:
        temp = visRange[0] + 24
        if temp > visRange[1]:
            v2 = [0, visRange[1]]
            v3 = [temp, 24]
        else:
            v2 = [0, 24]
    elif visRange[0] > 0 and visRange[1] > 24:
        temp = visRange[1] - 24
        if temp < visRange[0]:
            v1 = [0, temp]
            v2 = [visRange[0], 24]
        else:
            v2 = [0,24]
    else:
        v2 = [visRange[0],visRange[1]]
    
    visRanges = [v1,v2,v3]

    return visRanges

def inRange(visiblerange, ra): # Takes the visible range for a night and
    #the user's given RA and checks whether that RA is in the visible
    #range
    
    inrange = 0 # The given RA is outside the visible range (default)
    
    for subrange in visiblerange:
        if subrange != []:
            if ra in range(subrange[0],subrange[1]) or ra == subrange[1]:
                inrange = 1 # The given RA is in the visible range
    
    return inrange

#######################################################################  
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
########################### VISIBLE RANGES ############################
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
#######################################################################
"""
*** Goal ***
For the selected observatory, the user can input a date and get back 
what ranges of RA and DEC will be observable on that date.
"""

def getVisibleRange(observatory, time): # Takes the selected observatory 
    #and input date
    
    ### Check that the user provided all the correct inputs ###
    if observatory == "Select observatory":
        tk.messagebox.showerror("Error", "Please select an observatory.")
        return
        
    if time == "":
        tk.messagebox.showerror("Error", "Please enter a time.")
        return
        
    if len(time) < 10:
        tk.messagebox.showerror("Error", "Please format the date as yyyy-mm-dd.")
        return
        
    if time[4] != '-' or time[7] != '-':
        tk.messagebox.showerror("Error", "Please format the date as yyyy-mm-dd.")
        return
    
    ### Get and set observatory
    if observatory == 'CTIO (DECam, SOAR)':
        obsrvtry = astroplan.Observer.at_site("Cerro Tololo Interamerican Observatory",timezone='America/Santiago')
    elif observatory == 'VRO':
        obsrvtry = astroplan.Observer.at_site("Cerro Pachon",timezone='America/Santiago')
        tk.messagebox.showerror("Warning", "Horizon limits for VRO are currently unavailable. Using those for Blanco 4m")
    else:
        obsrvtry = astroplan.Observer.at_site("Kitt Peak National Observatory",timezone='US/Arizona')
        
    ### Determine sunrise and sunset in MST at the given observatory
    time = Time(time)
    
    sun_set = obsrvtry.sun_set_time(time, which="nearest") # In Julian
    sun_rise = obsrvtry.sun_rise_time(time, which="next") # In Julian
    
    sun_set = getDateTimefromJulian(sun_set) # UTC
    sun_rise = getDateTimefromJulian(sun_rise) #UTC
    
    yearSet=sun_set.split('-')[0]
    monthSet=sun_set.split('-')[1]
    daySet=sun_set.split('-')[2].split(' ')[0]
    hourSet=sun_set.split(' ')[1].split(':')[0]
    minuteSet=sun_set.split(' ')[1].split(':')[1]
    secSet=sun_set.split(' ')[1].split(':')[2]
    
    gmst = getGMST(yearSet, monthSet, 1, 0, 0, 0) # Setting day, hour
    #minute, and second to 1,0,0,0 to get GMST of the month
    
    mstSet = getMST(gmst, daySet, hourSet,minuteSet,secSet)
    
    yearRise=sun_rise.split('-')[0]
    monthRise=sun_rise.split('-')[1]
    dayRise=sun_rise.split('-')[2].split(' ')[0]
    hourRise=sun_rise.split(' ')[1].split(':')[0]
    minuteRise=sun_rise.split(' ')[1].split(':')[1]
    secRise=sun_rise.split(' ')[1].split(':')[2]
    
    gmst = getGMST(yearRise, monthRise, 1, 0, 0, 0) # Setting day, hour
    #minute, and second to 1,0,0,0 to get GMST of the month
    
    mstRise = getMST(gmst, dayRise, hourRise, minuteRise, secRise)
    
    ### Get the visible range for input night
    nightRange = getNightRange(mstSet,mstRise)
    
    tk.messagebox.showinfo("Visible Range for date", str(nightRange))
    
    return nightRange

#######################################################################  
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
########################### IS IT VISIBLE? ############################
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
#######################################################################
"""
*** Goal ***
The user inputs a set of RA and DEC coordinates, and the program returns 
whether or not the set is visible that day.
"""
def isVisible(observatory,ra,dec):
    
    ### Check that the user provided all the correct inputs ###
    if observatory == "Select observatory":
        tk.messagebox.showerror("Error", "Please select an observatory.")
        return
    
    ### Get and set observatory
    if observatory == 'CTIO (DECam, SOAR)':
        obsrvtry = astroplan.Observer.at_site("Cerro Tololo Interamerican Observatory",timezone='America/Santiago')
    elif observatory == 'VRO':
        obsrvtry = astroplan.Observer.at_site("Cerro Pachon",timezone='America/Santiago')
        tk.messagebox.showerror("Warning", "Horizon limits for VRO are currently unavailable. Using those for Blanco 4m")
    else:
        obsrvtry = astroplan.Observer.at_site("Kitt Peak National Observatory",timezone='US/Arizona')
       
    # Set default values for whether the user-given coordinates are visible
    raVis = 0 # RA not visible
    decVis = 1 # DEC visible
    
    if dec > 89 or dec < 37:
        decVis = 0 # In range
        
    # Get visible range for current date
    todayDate = Time(datetime.datetime.utcnow()) # Get the current date
    
    sun_set = obsrvtry.sun_set_time(todayDate, which="nearest") # In Julian
    sun_rise = obsrvtry.sun_rise_time(todayDate, which="next") # In Julian
    
    sun_set = getDateTimefromJulian(sun_set) # UTC
    sun_rise = getDateTimefromJulian(sun_rise) #UTC
    
    yearSet=sun_set.split('-')[0]
    monthSet=sun_set.split('-')[1]
    daySet=sun_set.split('-')[2].split(' ')[0]
    hourSet=sun_set.split(' ')[1].split(':')[0]
    minuteSet=sun_set.split(' ')[1].split(':')[1]
    secSet=sun_set.split(' ')[1].split(':')[2]
    
    gmst = getGMST(yearSet, monthSet, 1, 0, 0, 0) # Setting day, hour
    #minute, and second to 1,0,0,0 to get GMST of the month
    
    mstSet = getMST(gmst, daySet, hourSet,minuteSet,secSet)
    
    yearRise=sun_rise.split('-')[0]
    monthRise=sun_rise.split('-')[1]
    dayRise=sun_rise.split('-')[2].split(' ')[0]
    hourRise=sun_rise.split(' ')[1].split(':')[0]
    minuteRise=sun_rise.split(' ')[1].split(':')[1]
    secRise=sun_rise.split(' ')[1].split(':')[2]
    
    gmst = getGMST(yearRise, monthRise, 1, 0, 0, 0) # Setting day, hour
    #minute, and second to 1,0,0,0 to get GMST of the month
    
    mstRise = getMST(gmst, dayRise, hourRise, minuteRise, secRise)
        
    visibleRange = getNightRange(mstSet, mstRise)
    
    raVis = inRange(visibleRange, ra)
        
    if raVis + decVis == 2:
        yon = "Yes!"
    else:
        yon = "Nope."
    
    tk.messagebox.showinfo("Visible?", yon)
    
    return

#######################################################################  
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
############################### WIDGET ################################
#~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~8~#
#######################################################################

LARGE_FONT=("Verdant", 12)

# Observatories
observatories = [ "Select observatory", "CTIO (DECam, SOAR)", "VRO", "KNPO (DESI)" ] 

class widget(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        
        tk.Tk.wm_iconbitmap(self, "image-05Apr2020.ico")
        tk.Tk.wm_title(self, "Visible Sky")
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        
        self.frames = {}
        
        frame= MainPage(container, self)
        
        self.frames[MainPage] = frame
        
        frame.grid(row=0, column=0, sticky="nsew")
        
        self.show_frame(MainPage)
    
    def show_frame(self, cont):
        
        frame = self.frames[cont]
        frame.tkraise()
        
class MainPage(tk.Frame):
    
    def __init__(self, parent, controller):
        
        def visible():
            #print("This does nothing")
            getVisibleRange(variable.get(),user_time.get())
            
        def observableNow():
            
            if user_radec.get() == "" :
                tk.messagebox.showerror("Error", "Please enter an RA and DEC.")
                return
    
            ra,dec = user_radec.get().split(",")
            
            ra = float(ra)
            dec = float(dec)
            
            isVisible(variable.get(),ra,dec)
            
        
        tk.Frame.__init__(self, parent)
        
        label0 = tk.Label(self, text="Visible Sky", font=LARGE_FONT)
        label0.grid(pady=10, padx=10)
        
        variable = tk.StringVar(self)
        variable.set(observatories[0]) # default value
        
        w = ttk.OptionMenu(self, variable, *observatories) # Chosen Observatory
        w.grid(row=1,column=0)
        
        my_label1 = tk.Label(self, text = "Enter date in format yyyy-mm-dd ")
        my_label1.grid(row = 3, column = 0)
        user_time = tk.Entry(self) # User set time
        user_time.grid(row = 3, column = 1)
        
        my_button1 = ttk.Button(self, text = "Get visible ranges", command = visible)
        my_button1.grid(row = 3, column = 2)
        
        my_label2 = tk.Label(self, text = "Enter a RA and DEC in format RA(hr),DEC(deg)")
        my_label2.grid(row = 4, column = 0)
        user_radec = tk.Entry(self) # User set time
        user_radec.grid(row = 4, column = 1)
        
        my_button2 = ttk.Button(self, text = "Check if visible", command = observableNow)
        my_button2.grid(row = 4, column = 2)
        
app = widget()      
app.mainloop()
    