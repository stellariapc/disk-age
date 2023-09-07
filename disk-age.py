#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 00:48:23 2022

@author: mac
"""


#figsize
import numpy as np
import matplotlib.pyplot as plt
from xlrd import open_workbook
from scipy.spatial import ConvexHull
import math
from scipy.optimize import leastsq
import datetime
from scipy.optimize import minimize
from matplotlib.patches import Ellipse, Circle
from scipy.optimize import curve_fit
import csv
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import constants as const

from matplotlib.pyplot import MultipleLocator
from scipy.interpolate import splrep,splev,splint,sproot
import diptest
import scipy

def is_in_poly(p,poly):
    px,py = p
    is_in = False
    for i, corner in enumerate(poly):
        #print(i,corner)
        next_i = i + 1 if i + 1 < len(poly) else 0
        x1, y1 = corner
        x2, y2 = poly[next_i]
        if (x1 == px and y1 == py) or (x2 == px and y2 == py):
            is_in = True
            break
        if min(y1, y2) < py <= max(y1, y2):
            x = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
            if x == px:
                is_in = True
                break
            elif x > px:
                is_in = not is_in
    return is_in

def is_in_circle(p,circle):
    px,py = p
    center_x, center_y = circle[0]
    r = circle[1]
    if (px - center_x) ** 2 + (py - center_y) ** 2 <= r * r:
        return True
    else:
        return False

def is_in_ellipse(p,ellipse):
    px,py = p
    h,k,a,b,angle = ellipse
    #[srx,sry] 为沿p点顺时针旋转angle
    angle = np.radians(angle)
    srx = (px-h)*np.cos(angle) + (py-k)*np.sin(angle)+h
    sry = (py-k)*np.cos(angle) - (px-h)*np.sin(angle)+k
    if (srx-h)**2/a**2 + (sry-k)**2/b**2 <= 1:
        return True
    else:
        return False



def fun(x,k,b):
    return k*x+b

def get_XY_value_Taurus(list_magnitude,list_x,list_F0,A_AV,AV):
    #A2_AV,A3_AV,A4_AV,A5_AV,A8_AV = 0.112,0.058,0.023,0.022,0.020
    N_remove = []
    for i in range(len(list_magnitude)):
        if list_magnitude[i] == '':
            N_remove.append(i)
    N_remove.reverse()
    #print(list_magnitude,N_remove)
    for i in N_remove:
        list_magnitude.pop(i)
        list_x.pop(i)
        list_F0.pop(i)
        A_AV.pop(i)
    list_magnitude_revised = [float(list_magnitude[j])-A_AV[j]*AV for j in range(len(list_magnitude))]
    Y = [list_F0[i]-list_magnitude_revised[i]/2.5 - list_x[i]  for i in range(len(list_magnitude_revised))]
    return list_x,Y

def get_XY_value_Orion(list_magnitude,list_x,list_F0,A_AK,AK):
    #A2_AV,A3_AV,A4_AV,A5_AV,A8_AV = 0.112,0.058,0.023,0.022,0.020
    N_remove = []
    for i in range(len(list_magnitude)):
        if list_magnitude[i] == '':
            N_remove.append(i)
    N_remove.reverse()
    for i in N_remove:
        list_magnitude.pop(i)
        list_x.pop(i)
        list_F0.pop(i)
        A_AK.pop(i)
    list_magnitude_revised = [float(list_magnitude[j])-A_AK[j]*AK for j in range(len(list_magnitude))]
    Y = [list_F0[i]-list_magnitude_revised[i]/2.5 - list_x[i] for i in range(len(list_magnitude_revised))]
    return list_x,Y

def get_XY_value_15(list_F,list_x,A_AV,AV):
    N_remove = []
    for i in range(len(list_F)):
        if list_F[i] == '':
            N_remove.append(i)
    N_remove.reverse()
    for i in N_remove:
        list_F.pop(i)
        list_x.pop(i)
        A_AV.pop(i)
    Y = [np.log10(float(list_F[j]) * ((100**(1/5))**(A_AV[j]*AV)) ) - list_x[j] for j in range(len(list_F))]#list_F单位是mJY/F_Hz,需要转换成F_lamda
    return list_x,Y   


def Get_coor_alpha_with_24_no_revise(Protocluster,num,Taurus_known,Taurus_new,Orion,PCs_18):
    coor,alpha_PC = [],[]
    DATA = []
    DATA_2MASS = []
    if Protocluster[num]['name'] == 'Taurus':
        csvFile = open('Taurus_known 2MASS.csv', "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA_2MASS.append(row)
        csvFile = open('Taurus_new 2MASS.csv', "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA_2MASS.append(row)
        
        csvFile = open(Taurus_known, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        csvFile = open(Taurus_new, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            x,y = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if x >= Protocluster[num]['field'][0] and x <= Protocluster[num]['field'][1] \
            and y >= Protocluster[num]['field'][2] and y <= Protocluster[num]['field'][3]:
                K = ''
                for j in range(len(DATA_2MASS)):
                    if [x,y] == [float(DATA_2MASS[j]['_RAJ2000']),float(DATA_2MASS[j]['_DEJ2000'])] :
                        J,H,K = DATA_2MASS[j]['Jmag'],DATA_2MASS[j]['Hmag'],DATA_2MASS[j]['Kmag']
                        break
                    
                lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['[3.6]'],DATA[i]['[4.5]'],DATA[i]['[5.8]'],DATA[i]['[8.0]'],DATA[i]['[24]']
                
                list_magnitude = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]
                list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]
                list_F0 = [np.log10(666.7),
                         np.log10(280.9),np.log10(179.7),
                         np.log10(115.0),np.log10(64.9),
                         np.log10(7.17) ] #F_0
            
                if DATA[i]['AV'] == '':
                    AV = Protocluster[num]['AV_average']
                else:
                    AV = float(DATA[i]['AV'])
                A_AK = [1.0,
                         0.5916522648641053,
                         0.4854745133467907,
                         0.39437870467923686,
                         0.450759293776694,
                         0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                A_AV = [0.11330520669381368,
                         0.06703728216129044,
                         0.055006790079336725,
                         0.04468516064931943,
                         0.051073374950525795,
                         0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction

                X,Y = get_XY_value_Taurus(list_magnitude,list_x,list_F0,A_AV,AV)
                
                if len(X) >= 2 and lamda_24 != '':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    alpha_PC.append(alpha)
                    coor.append([x,y])
    elif Protocluster[num]['name'] in ['Orion A','Orion B']:
        csvFile = open(Orion, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            x,y = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if x >= Protocluster[num]['field'][0] and x <= Protocluster[num]['field'][1] \
            and y >= Protocluster[num]['field'][2] and y <= Protocluster[num]['field'][3]:
                J,H,K = DATA[i]['Jmag'],DATA[i]['Hmag'],DATA[i]['Kmag']
                lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['[3.6]'],DATA[i]['[4.5]'],DATA[i]['[5.8]'],DATA[i]['[8.0]'],DATA[i]['[24]']
                
                list_magnitude = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]
                list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]
                list_F0 = [np.log10(666.7),
                         np.log10(280.9),np.log10(179.7),
                         np.log10(115.0),np.log10(64.9),
                         np.log10(7.17) ] #F_0
                if DATA[i]['AKs'] == '':
                    AK = Protocluster[num]['AK_average']
                else:
                    AK = float(DATA[i]['AKs'])
                A_AK = [1.0,
                         0.5916522648641053,
                         0.4854745133467907,
                         0.39437870467923686,
                         0.450759293776694,
                         0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                A_AV = [0.11330520669381368,
                         0.06703728216129044,
                         0.055006790079336725,
                         0.04468516064931943,
                         0.051073374950525795,
                         0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction
                X,Y = get_XY_value_Orion(list_magnitude,list_x,list_F0,A_AK,AK)
                if len(X) >= 2 and lamda_24 != '':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    alpha_PC.append(alpha)
                    coor.append([x,y])
                
    else:
        csvFile = open(PCs_18, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            if DATA[i]['Cloud'] == Protocluster[num]['name']:
                RA, DEC = float(DATA[i]['_RAJ2000']), float(DATA[i]['_DEJ2000'])
                if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] and \
                   DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                     
                    
                    J,H,K = DATA[i]['F1.25'],DATA[i]['F1.65'],DATA[i]['F2.17']
                    lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['F3.6'],DATA[i]['F4.5'],DATA[i]['F5.8'],DATA[i]['F8.0'],DATA[i]['F24']
                    lamda_12,lamda_14,lamda_22 = DATA[i]['F12'],DATA[i]['F14'],DATA[i]['F22']
                    lamda_3_4,lamda_5,lamda_6_7 = DATA[i]['F3.4'],DATA[i]['F5.0'],DATA[i]['F6.7']
                    
                    list_F = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]#
                    list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]#
                    AV = float(DATA[i]['Av'])
                    A_AK = [1.0,
                             0.5916522648641053,
                             0.4854745133467907,
                             0.39437870467923686,
                             0.450759293776694,
                             0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                    A_AV = [0.11330520669381368,
                             0.06703728216129044,
                             0.055006790079336725,
                             0.04468516064931943,
                             0.051073374950525795,
                             0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction

                    X,Y = get_XY_value_15(list_F,list_x,A_AV,AV)
                    if len(X) >= 2 and DATA[i]['AGB?'] == 'N' \
                        and lamda_24 != '':
                        k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                        b0 = Y[0] - k0 * X[0]
                        popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                        alpha = popt[0]
                        alpha_PC.append(alpha)
                        coor.append([RA, DEC])
    return coor,alpha_PC


def Get_coor_alpha_with_24(Protocluster,num,Taurus_known,Taurus_new,Orion,PCs_18):
    coor,alpha_PC = [],[]
    DATA = []
    DATA_2MASS = []
    if Protocluster[num]['name'] == 'Taurus':
        csvFile = open('Taurus_known 2MASS.csv', "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA_2MASS.append(row)
        csvFile = open('Taurus_new 2MASS.csv', "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA_2MASS.append(row)
        
        csvFile = open(Taurus_known, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        csvFile = open(Taurus_new, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            x,y = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if x >= Protocluster[num]['field'][0] and x <= Protocluster[num]['field'][1] \
            and y >= Protocluster[num]['field'][2] and y <= Protocluster[num]['field'][3]:
                K = ''
                for j in range(len(DATA_2MASS)):
                    if [x,y] == [float(DATA_2MASS[j]['_RAJ2000']),float(DATA_2MASS[j]['_DEJ2000'])] :
                        J,H,K = DATA_2MASS[j]['Jmag'],DATA_2MASS[j]['Hmag'],DATA_2MASS[j]['Kmag']
                        break
                    
                lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['[3.6]'],DATA[i]['[4.5]'],DATA[i]['[5.8]'],DATA[i]['[8.0]'],DATA[i]['[24]']
                
                list_magnitude = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]
                list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]
                list_F0 = [np.log10(666.7),
                         np.log10(280.9),np.log10(179.7),
                         np.log10(115.0),np.log10(64.9),
                         np.log10(7.17) ] #F_0
            
                if DATA[i]['AV'] == '':
                    AV = Protocluster[num]['AV_average']
                else:
                    AV = float(DATA[i]['AV'])
                A_AK = [1.0,
                         0.5916522648641053,
                         0.4854745133467907,
                         0.39437870467923686,
                         0.450759293776694,
                         0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                A_AV = [0.11330520669381368,
                         0.06703728216129044,
                         0.055006790079336725,
                         0.04468516064931943,
                         0.051073374950525795,
                         0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction

                X,Y = get_XY_value_Taurus(list_magnitude,list_x,list_F0,A_AV,AV)
                
                if len(X) >= 2 and lamda_24 != '':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    alpha_PC.append(alpha)
                    coor.append([x,y])
    elif Protocluster[num]['name'] in ['Orion A','Orion B']:
        csvFile = open(Orion, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            x,y = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if x >= Protocluster[num]['field'][0] and x <= Protocluster[num]['field'][1] \
            and y >= Protocluster[num]['field'][2] and y <= Protocluster[num]['field'][3]:
                J,H,K = DATA[i]['Jmag'],DATA[i]['Hmag'],DATA[i]['Kmag']
                lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['[3.6]'],DATA[i]['[4.5]'],DATA[i]['[5.8]'],DATA[i]['[8.0]'],DATA[i]['[24]']
                
                list_magnitude = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]
                list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]
                list_F0 = [np.log10(666.7),
                         np.log10(280.9),np.log10(179.7),
                         np.log10(115.0),np.log10(64.9),
                         np.log10(7.17) ] #F_0
                if DATA[i]['AKs'] == '':
                    AK = Protocluster[num]['AK_average']
                else:
                    AK = float(DATA[i]['AKs'])
                A_AK = [1.0,
                         0.5916522648641053,
                         0.4854745133467907,
                         0.39437870467923686,
                         0.450759293776694,
                         0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                A_AV = [0.11330520669381368,
                         0.06703728216129044,
                         0.055006790079336725,
                         0.04468516064931943,
                         0.051073374950525795,
                         0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction
                X,Y = get_XY_value_Orion(list_magnitude,list_x,list_F0,A_AK,AK)
                if len(X) >= 2 and lamda_24 != '':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    alpha_PC.append(alpha)
                    coor.append([x,y])
                
                if len(X) >= 2 and lamda_24 == '' and lamda_8 != '' and \
                    Protocluster[num]['savename'] != 'Orion B-N':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    if Protocluster[num]['savename'] == 'Orion A':
                        k,b = 0.5861471324025159, -0.06808291758577195
                    if Protocluster[num]['savename'] == 'Orion B-M':
                        k,b = 0.6630901796382047, 0.007461135521107512
                    if Protocluster[num]['savename'] == 'Orion B-S':
                        k,b = 0.613707133224276 , -0.09051645461742297
                    alpha_revised = k * alpha + b
                    alpha_PC.append(alpha_revised)
                    coor.append([x,y])
                
                if len(X) >= 2 and lamda_24 == '' and lamda_8 == '' and lamda_5_8 != '' and \
                    Protocluster[num]['savename'] != 'Orion B-N':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    if Protocluster[num]['savename'] == 'Orion A':
                        k,b = 0.4791595379570749, -0.1328556839959108
                    if Protocluster[num]['savename'] == 'Orion B-M':
                        k,b = 0.5521702053055223, -0.04322091433810932
                    if Protocluster[num]['savename'] == 'Orion B-S':
                        k,b = 0.508129372348308, -0.12814708605229747
                    alpha_revised = k * alpha + b
                    alpha_PC.append(alpha_revised)
                    coor.append([x,y])
                
    else:
        csvFile = open(PCs_18, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            if DATA[i]['Cloud'] == Protocluster[num]['name']:
                RA, DEC = float(DATA[i]['_RAJ2000']), float(DATA[i]['_DEJ2000'])
                if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] and \
                   DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                     
                    
                    J,H,K = DATA[i]['F1.25'],DATA[i]['F1.65'],DATA[i]['F2.17']
                    lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['F3.6'],DATA[i]['F4.5'],DATA[i]['F5.8'],DATA[i]['F8.0'],DATA[i]['F24']
                    lamda_12,lamda_14,lamda_22 = DATA[i]['F12'],DATA[i]['F14'],DATA[i]['F22']
                    lamda_3_4,lamda_5,lamda_6_7 = DATA[i]['F3.4'],DATA[i]['F5.0'],DATA[i]['F6.7']
                    
                    list_F = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]#
                    list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]#
                    AV = float(DATA[i]['Av'])
                    A_AK = [1.0,
                             0.5916522648641053,
                             0.4854745133467907,
                             0.39437870467923686,
                             0.450759293776694,
                             0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                    A_AV = [0.11330520669381368,
                             0.06703728216129044,
                             0.055006790079336725,
                             0.04468516064931943,
                             0.051073374950525795,
                             0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction

                    X,Y = get_XY_value_15(list_F,list_x,A_AV,AV)
                    if len(X) >= 2 and DATA[i]['AGB?'] == 'N' \
                        and lamda_24 != '':
                        k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                        b0 = Y[0] - k0 * X[0]
                        popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                        alpha = popt[0]
                        alpha_PC.append(alpha)
                        coor.append([RA, DEC])
    return coor,alpha_PC


#交叉认证可能会在给定范围中对一个源证认出多个源
def get_distance(Protocluster,num):
    coor,coor_gaia = [],[]
    dis,dis_low,dis_up = [],[],[]
    id_source = []
    DATA = []
    if Protocluster[num]['name'] not in ['Orion A','Orion B','Taurus']:
        csvFile = open("18 in Gaia EDR3 distance 3 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            if DATA[i]['Cloud'] == Protocluster[num]['name']:
                RA, DEC = float(DATA[i]['_RAJ2000']), float(DATA[i]['_DEJ2000'])
                if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] and \
                   DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                    if DATA[i]['AGB?']=='N':
                        if DATA[i]['rgeo'] != '':
                            coor.append([RA, DEC])
                            id_source.append(DATA[i]['Source'])
                            coor_gaia.append([float(DATA[i]['RA_ICRS']), float(DATA[i]['DE_ICRS'])])
                            dis.append(float(DATA[i]['rgeo']))
                            dis_low.append(float(DATA[i]['b_rgeo']))
                            dis_up.append(float(DATA[i]['B_rgeo']))
    if Protocluster[num]['name'] in ['Orion A','Orion B']:
        csvFile = open("Orion in Gaia EDR3 distance 3 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            RA, DEC = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] \
            and DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                if DATA[i]['rgeo'] != '':
                    coor.append([RA, DEC])
                    id_source.append(DATA[i]['Source'])
                    coor_gaia.append([float(DATA[i]['RA_ICRS']), float(DATA[i]['DE_ICRS'])])
                    dis.append(float(DATA[i]['rgeo']))
                    dis_low.append(float(DATA[i]['b_rgeo']))
                    dis_up.append(float(DATA[i]['B_rgeo']))
    if Protocluster[num]['name'] in ['Taurus']:
        csvFile = open("Taurus_known in Gaia EDR3 distance 3 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        csvFile = open("Taurus_new in Gaia EDR3 distance 3 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            RA, DEC = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] \
            and DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                if DATA[i]['rgeo'] != '':
                    coor.append([RA, DEC])
                    id_source.append(DATA[i]['Source'])
                    coor_gaia.append([float(DATA[i]['RA_ICRS']), float(DATA[i]['DE_ICRS'])])
                    dis.append(float(DATA[i]['rgeo']))
                    dis_low.append(float(DATA[i]['b_rgeo']))
                    dis_up.append(float(DATA[i]['B_rgeo']))
    return coor,coor_gaia,dis,dis_low,dis_up,id_source


def get_pm_parallax(Protocluster,num):
    coor_spitzer = []
    pm_ra,pm_ra_error,pm_dec,pm_dec_error = [],[],[],[]
    parallax,parallax_error = [],[]
    id_source = []
    DATA = []
    if Protocluster[num]['name'] not in ['Orion A','Orion B','Taurus']:
        csvFile = open("18 in Gaia EDR3 2 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            if DATA[i]['Cloud'] == Protocluster[num]['name']:
                RA, DEC = float(DATA[i]['_RAJ2000']), float(DATA[i]['_DEJ2000'])
                if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] and \
                   DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                    if DATA[i]['AGB?']=='N':
                        if DATA[i]['pmra'] != '':
                            if i == 0:
                                coor_spitzer.append([RA, DEC])
                                id_source.append(DATA[i]['source_id'])
                                pm_ra.append(float(DATA[i]['pmra']))
                                pm_ra_error.append(float(DATA[i]['pmra_error']))
                                pm_dec.append(float(DATA[i]['pmdec']))
                                pm_dec_error.append(float(DATA[i]['pmdec_error']))
                                parallax.append(float(DATA[i]['parallax']))
                                parallax_error.append(float(DATA[i]['parallax_error']))
                            if i > 0:
                                if [RA, DEC] != [float(DATA[i-1]['_RAJ2000']), float(DATA[i-1]['_DEJ2000'])]:
                                    coor_spitzer.append([RA, DEC])
                                    id_source.append(DATA[i]['source_id'])
                                    pm_ra.append(float(DATA[i]['pmra']))
                                    pm_ra_error.append(float(DATA[i]['pmra_error']))
                                    pm_dec.append(float(DATA[i]['pmdec']))
                                    pm_dec_error.append(float(DATA[i]['pmdec_error']))
                                    parallax.append(float(DATA[i]['parallax']))
                                    parallax_error.append(float(DATA[i]['parallax_error']))
    if Protocluster[num]['name'] in ['Orion A','Orion B']:
        csvFile = open("Orion in Gaia EDR3 2 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            RA, DEC = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] \
            and DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                if DATA[i]['pmra'] != '':
                    if i == 0:
                        coor_spitzer.append([RA, DEC])
                        id_source.append(DATA[i]['source_id'])
                        pm_ra.append(float(DATA[i]['pmra']))
                        pm_ra_error.append(float(DATA[i]['pmra_error']))
                        pm_dec.append(float(DATA[i]['pmdec']))
                        pm_dec_error.append(float(DATA[i]['pmdec_error']))
                        parallax.append(float(DATA[i]['parallax']))
                        parallax_error.append(float(DATA[i]['parallax_error']))
                    if i > 0:
                        if [RA, DEC] != [float(DATA[i-1]['_RAJ2000']), float(DATA[i-1]['_DEJ2000'])]:
                            coor_spitzer.append([RA, DEC])
                            id_source.append(DATA[i]['source_id'])
                            pm_ra.append(float(DATA[i]['pmra']))
                            pm_ra_error.append(float(DATA[i]['pmra_error']))
                            pm_dec.append(float(DATA[i]['pmdec']))
                            pm_dec_error.append(float(DATA[i]['pmdec_error']))
                            parallax.append(float(DATA[i]['parallax']))
                            parallax_error.append(float(DATA[i]['parallax_error']))
    if Protocluster[num]['name'] in ['Taurus']:
        csvFile = open("Taurus_known in Gaia EDR3 2 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        csvFile = open("Taurus_new in Gaia EDR3 2 arcsec.csv", "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            RA, DEC = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] \
            and DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                if DATA[i]['pmra'] != '':
                    if i == 0:
                        coor_spitzer.append([RA, DEC])
                        id_source.append(DATA[i]['source_id'])
                        pm_ra.append(float(DATA[i]['pmra']))
                        pm_ra_error.append(float(DATA[i]['pmra_error']))
                        pm_dec.append(float(DATA[i]['pmdec']))
                        pm_dec_error.append(float(DATA[i]['pmdec_error']))
                        parallax.append(float(DATA[i]['parallax']))
                        parallax_error.append(float(DATA[i]['parallax_error']))
                    if i > 0:
                        if [RA, DEC] != [float(DATA[i-1]['_RAJ2000']), float(DATA[i-1]['_DEJ2000'])]:
                            coor_spitzer.append([RA, DEC])
                            id_source.append(DATA[i]['source_id'])
                            pm_ra.append(float(DATA[i]['pmra']))
                            pm_ra_error.append(float(DATA[i]['pmra_error']))
                            pm_dec.append(float(DATA[i]['pmdec']))
                            pm_dec_error.append(float(DATA[i]['pmdec_error']))
                            parallax.append(float(DATA[i]['parallax']))
                            parallax_error.append(float(DATA[i]['parallax_error']))
    return coor_spitzer,pm_ra,pm_ra_error,pm_dec,pm_dec_error,parallax,parallax_error,id_source



def Get_coor_alpha_with_24_flux(Protocluster,num,Taurus_known,Taurus_new,Orion,PCs_18):
    FLUX_x = []
    FLUX_y = []
    coor,alpha_PC = [],[]
    DATA = []
    DATA_2MASS = []
    if Protocluster[num]['name'] == 'Taurus':
        csvFile = open('Taurus_known 2MASS.csv', "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA_2MASS.append(row)
        csvFile = open('Taurus_new 2MASS.csv', "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA_2MASS.append(row)
        
        csvFile = open(Taurus_known, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        csvFile = open(Taurus_new, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            x,y = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if x >= Protocluster[num]['field'][0] and x <= Protocluster[num]['field'][1] \
            and y >= Protocluster[num]['field'][2] and y <= Protocluster[num]['field'][3]:
                K = ''
                for j in range(len(DATA_2MASS)):
                    if [x,y] == [float(DATA_2MASS[j]['_RAJ2000']),float(DATA_2MASS[j]['_DEJ2000'])] :
                        J,H,K = DATA_2MASS[j]['Jmag'],DATA_2MASS[j]['Hmag'],DATA_2MASS[j]['Kmag']
                        break
                    
                lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['[3.6]'],DATA[i]['[4.5]'],DATA[i]['[5.8]'],DATA[i]['[8.0]'],DATA[i]['[24]']
                
                list_magnitude = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]
                list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]
                list_F0 = [np.log10(666.7),
                         np.log10(280.9),np.log10(179.7),
                         np.log10(115.0),np.log10(64.9),
                         np.log10(7.17) ] #F_0
            
                if DATA[i]['AV'] == '':
                    AV = Protocluster[num]['AV_average']
                else:
                    AV = float(DATA[i]['AV'])
                A_AK = [1.0,
                         0.5916522648641053,
                         0.4854745133467907,
                         0.39437870467923686,
                         0.450759293776694,
                         0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                A_AV = [0.11330520669381368,
                         0.06703728216129044,
                         0.055006790079336725,
                         0.04468516064931943,
                         0.051073374950525795,
                         0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction

                X,Y = get_XY_value_Taurus(list_magnitude,list_x,list_F0,A_AV,AV)
                
                if K != '' and lamda_24 != '':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    alpha_PC.append(alpha)
                    coor.append([x,y])
                    FLUX_x.append(X)
                    FLUX_y.append(Y)
    elif Protocluster[num]['name'] in ['Orion A','Orion B']:
        csvFile = open(Orion, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            x,y = float(DATA[i]['_RAJ2000']),float(DATA[i]['_DEJ2000']) 
            if x >= Protocluster[num]['field'][0] and x <= Protocluster[num]['field'][1] \
            and y >= Protocluster[num]['field'][2] and y <= Protocluster[num]['field'][3]:
                J,H,K = DATA[i]['Jmag'],DATA[i]['Hmag'],DATA[i]['Kmag']
                lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['[3.6]'],DATA[i]['[4.5]'],DATA[i]['[5.8]'],DATA[i]['[8.0]'],DATA[i]['[24]']
                
                list_magnitude = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]
                list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]
                list_F0 = [np.log10(666.7),
                         np.log10(280.9),np.log10(179.7),
                         np.log10(115.0),np.log10(64.9),
                         np.log10(7.17) ] #F_0
                if DATA[i]['AKs'] == '':
                    AK = Protocluster[num]['AK_average']
                else:
                    AK = float(DATA[i]['AKs'])
                A_AK = [1.0,
                         0.5916522648641053,
                         0.4854745133467907,
                         0.39437870467923686,
                         0.450759293776694,
                         0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                A_AV = [0.11330520669381368,
                         0.06703728216129044,
                         0.055006790079336725,
                         0.04468516064931943,
                         0.051073374950525795,
                         0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction
                X,Y = get_XY_value_Orion(list_magnitude,list_x,list_F0,A_AK,AK)
                if K != '' and lamda_24 != '':
                    k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                    b0 = Y[0] - k0 * X[0]
                    popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                    alpha = popt[0]
                    alpha_PC.append(alpha)
                    coor.append([x,y])
                    FLUX_x.append(X)
                    FLUX_y.append(Y)
                
    else:
        csvFile = open(PCs_18, "r")
        dict_reader = csv.DictReader(csvFile)
        for row in dict_reader:
            DATA.append(row)
        for i in range(len(DATA)):
            if DATA[i]['Cloud'] == Protocluster[num]['name']:
                RA, DEC = float(DATA[i]['_RAJ2000']), float(DATA[i]['_DEJ2000'])
                if RA >= Protocluster[num]['field'][0] and RA <= Protocluster[num]['field'][1] and \
                   DEC >= Protocluster[num]['field'][2] and DEC <= Protocluster[num]['field'][3]:
                    
                    J,H,K = DATA[i]['F1.25'],DATA[i]['F1.65'],DATA[i]['F2.17']
                    lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24 = DATA[i]['F3.6'],DATA[i]['F4.5'],DATA[i]['F5.8'],DATA[i]['F8.0'],DATA[i]['F24']
                    lamda_12,lamda_14,lamda_22 = DATA[i]['F12'],DATA[i]['F14'],DATA[i]['F22']
                    lamda_3_4,lamda_5,lamda_6_7 = DATA[i]['F3.4'],DATA[i]['F5.0'],DATA[i]['F6.7']
                    
                    list_F = [K,lamda_3_6,lamda_4_5,lamda_5_8,lamda_8,lamda_24]#
                    list_x = [np.log10(2.159),np.log10(3.6),np.log10(4.5),np.log10(5.8),np.log10(8),np.log10(23.68)]#
                    AV = float(DATA[i]['Av'])
                    A_AK = [1.0,
                             0.5916522648641053,
                             0.4854745133467907,
                             0.39437870467923686,
                             0.450759293776694,
                             0.23709433802459026] #only contain 2.17,3.6,4.5,5.8,8.0um,at last, 24um, dont care extinction
                    A_AV = [0.11330520669381368,
                             0.06703728216129044,
                             0.055006790079336725,
                             0.04468516064931943,
                             0.051073374950525795,
                             0.026864022975809126] #only contain 2.17,3.6,4.5,5.8,8.0,,at last,24um, dont care extinction

                    X,Y = get_XY_value_15(list_F,list_x,A_AV,AV)
                    if len(X) >= 2 and DATA[i]['AGB?'] == 'N' \
                        and lamda_24 != '' and K != '':
                        k0 = ( Y[len(X)-1] - Y[0] )/( X[len(X)-1] - X[0] )
                        b0 = Y[0] - k0 * X[0]
                        popt,pcov = curve_fit(fun, X, Y, p0=[k0,b0])
                        alpha = popt[0]
                        alpha_PC.append(alpha)
                        coor.append([RA, DEC])
                        FLUX_x.append(X)
                        FLUX_y.append(Y)
    return coor,alpha_PC,FLUX_x,FLUX_y




#[274,  280,     -5,     -1], 

Protocluster = [
                {'name':'Aquila',             'savename':'Serpens NE',        
                 'field':[278.75,281,     -1,      1],
                 'distance_select':[400, 600] },
                
                {'name':'Aquila',             'savename':'Aquila-Main',        
                 'field':[275.5,  278.6,     -4.3,     -1.2], 
                 'distance_select':[400, 600] },
                 
                {'name':'Ophiuchus',              'savename':'Ophiuchus',     
                 'field':[245,   252,     -27,   -21],
                 'distance_select':[100,200] },
                
                {'name':'Serpens',                'savename':'Serpens-Main',       
                 'field':[276,  278,    -1,     2],   
                 'distance_select':[300,600] },
                
                {'name':'Chamaeleon I',           'savename':'Chamaeleon I',  
                 'field':[163,  170,    -79,    -75],
                 'distance_select':[150,300]},
                
                {'name':'Perseus',                'savename':'Perseus-E',       
                 'field':[54,   58,     29,       34],
                 'distance_select':[200,400] },
                
                {'name':'Taurus',                 'savename':'Taurus',        
                 'field':[62,  75,  20, 32], 
                 'AV_average':3.808,
                 'distance_select':[100,200] },
                
                
                
                {'name':'IC5146',                 'savename':'IC 5146-E',        
                 'field':[327.5,329,     46,     49],
                 'distance_select':[600,900] },
                
                
                
                
                {'name':'Orion A',                'savename':'Orion A',        
                 'field':[82.6,   88,     -11,   -3.5], 
                 'AK_average':0.394,
                 'distance_select':[300,600] },
                
                {'name':'Orion B',                'savename':'Orion B-S',        
                 'field':[85,   86,     -3.5,   2.1], 
                 'AK_average':0.871,
                 'distance_select':[300,600] },
                
                {'name':'Orion B',                'savename':'Orion B-M',        
                 'field':[86,   88,     -0.7,   2.1], 
                 'AK_average':0.612,
                 'distance_select':[300,600] },
                
                
                
                
                {'name':'Auriga/CMC',             'savename':'AurigaCMC',     
                 'field':[60 ,  70,      34,     42],
                 'distance_select':[400, 600] },
                
                {'name':'Perseus',                'savename':'Perseus-W',       
                 'field':[50,   54,     29,       34],
                 'distance_select':[200,400] },
                
                
                
                {'name':'Orion B',                'savename':'Orion B-N',        
                 'field':[88,   90,     1,   2.1], 
                 'AK_average':0.113,
                 'distance_select':[200,400] },
                
                {'name':'Cepheus',               'savename':'Cepheus-SW',       
                 'field':[312,  320,     60,     70],
                 'distance_select':[300,400] },
                {'name':'Cepheus',               'savename':'Cepheus-NW',       
                 'field':[312,  320,     70,     80],
                 'distance_select':[300,400] },
                {'name':'Cepheus',               'savename':'Cepheus-NE',       
                 'field':[335,  345,     70,     80],
                 'distance_select':[300,400] },
                
                {'name':'Chamaeleon II',         'savename':'Chamaeleon II', 
                 'field':[191,  198,    -80,    -75],
                 'distance_select':[150,300] },
                
                {'name':'Corona Australis',      'savename':'Corona Australis',
                 'field':[284,  287,    -38,    -36],
                 'distance_select':[100,200] },
                {'name':'IC5146',                'savename':'IC 5146-W',        
                 'field':[325,  327.5,   46,     49],
                 'distance_select':[600,900] },
                
                {'name':'Lupus I',               'savename':'Lupus I',       
                 'field':[0,   360,     -90,     90],
                 'distance_select':[100,200] },
                {'name':'Lupus III',             'savename':'Lupus III',     
                 'field':[241,  244,     -40,     -38],
                 'distance_select':[100,200] },
                {'name':'Lupus IV',              'savename':'Lupus IV',      
                 'field':[0,   360,     -90,     90],
                 'distance_select':[100,200] },
                {'name':'Lupus V',               'savename':'Lupus V',       
                 'field':[244,  247,     -39,     -36],
                 'distance_select':[100,200] },
                {'name':'Lupus VI',              'savename':'Lupus VI',      
                 'field':[244,  248,     -42,     -39],
                 'distance_select':[100,200] },
                {'name':'Musca',                 'savename':'Musca',         
                 'field':[0,   360,     -90,     90],
                 'distance_select':[200,600] },
                
                {'name':'Ophiuchus North',       'savename':'Ophiuchus North',    
                 'field':[0,   360,     -90,     90], 
                 'distance_select':[100,200] },
                
                
                {'name':'Chamaeleon III',        'savename':'Chamaeleon III',
                 'field':[0,   360,     -90,     90],
                 'distance_select':[150,300]},

                ]


Protocluster_cared = [
                     {'name':'Ophiuchus',      'savename':'Ophiuchus', 
                      'distance_error':[138,6],
                      'innerarea':[[[252,-27],[245,-27],[245,-21],[252,-21],[252,-27]],
                                   [[246.75,-24.5],0.35],],
                      },


                      
                     {'name':'Perseus',    'savename':'Perseus-E', 
                      'distance_error':[308,21],

                      'innerarea':[[[57,31],[54.5,31],[54.5,33],[57,33],[57,31]],
                                   [[56.15,32.21,0.15,0.12,80]],
                                   [[56.02,32.04,0.18,0.07,-20]],
                                   
                                   [[55.67,31.79,0.6,0.3,-16]],
                                   ],
                      },
                     

                      
                     {'name':'Orion B',      'savename':'Orion B-M', 
                     'distance_error':[403,40],
                     'innerarea':[[[87.2,-0.6],[86.2,-0.6],[86.2,0.8],[87.2,0.8],[87.2,-0.6]],
                                  [[86.80,0.38],0.14],
                                  [[86.65,0.15],0.12],
                                  [[86.6,-0.12,0.17,0.1,60]],
                                  ],
                     
                     }, 
                      
                      
                      
                    {'name':'Orion B',        'savename':'Orion B-S',  
                     'distance_error':[391,45],

                     'innerarea':[[[86.25,-2.8],[85,-2.8],[85,-1],[86.25,-1],[86.25,-2.8]],
                                  [[85.43,-1.85],0.25],
                                  [[85.4,-2.3],0.15],
                                  ],
                     
                     },
                      
                     
                    {'name':'Perseus',       'savename':'Perseus-W',  
                     'distance_error':[289,21],
                     
                     'innerarea':[[[54,30],[51,30],[51,32],[54,32],[54,30]],
                                  
                                  [[52.25,31.3,0.3,0.2,60]],
                                  [[53.0,30.85,0.8,0.4,45]],],
                     
                     },
                     
                    {'name':'Serpens',         'savename':'Serpens-Main',  
                      'distance_error':[433,39],
                      
                      'innerarea':[[[278,-0.5],[276.75,-0.5],[276.75,1.5],[278,1.5],[278,-0.5]],
                                   [[277.5,1.25],0.1],
                                   [[277.45,0.62,0.2,0.13,80]],
                                   [[277.25,0.5],0.1],
                                   [[277.25,0.3],0.1],],
                      },
                     
                     
                     {'name':'Taurus',           'savename':'Taurus',      
                      'distance_error':[137,14],
                       'innerarea':[[[73,21],[63,21],[63,30],[73,30],[73,21]],
                                    
                                    #[[69.3,26.3],[68,25],[61,29.5],[64.3,29.5],[69.3,26.3]],
                                    [[68,25],[62,28],[64,30],[70,27],[68,25]],
                                    [[69.3,26.3],[69.3,25],[71,25],[71,26.3],[69.3,26.3]],
                                    #[[66,23.3],[66,25],[71,25],[71,23.3],[66,23.3]],
                                    
                                    [[64,24.7],[64,25.8],[71,24],[71,23.3],[67,23.4],[64,24.7]],
                                    
                                    [[67.5,22],[67.5,23.3],[70,23.3],[70,22],[67.5,22]],
                                    
                                    ],
                       },
                     
                     {'name':'IC5146',            'savename':'IC 5146-E',    
                      'distance_error':[754,48],
                        'innerarea':[[[328.7,47.1],[328,47.1],[328,47.5],[328.7,47.5],[328.7,47.1]],
                                     [[328.40,47.27],0.11],
                                     [[328.17,47.2],0.06],
                                     ],
                       },
                     
                     
                     
                     
                     
                     {'name':'Chamaeleon I',      'savename':'Chamaeleon I', 
                        'innerarea':[[[163,-78],[169,-78],[169,-76],[163,-76],[163,-78]],
                                     [[166,-76],[168.5,-76],[168.5,-77],[166,-77],[166,-76]],
                                     [[166,-78],[168.5,-78],[168.5,-77],[166,-77],[166,-78]],
                                     ],
                        
                        },
                     
                    {'name':'Aquila',      'savename':'Aquila-Main',
                     'innerarea':[
                                 [[275.5,-4.3],[275.5,-1.2],[278.6,-1.2],[278.6,-4.3],[275.5,-4.3]],
                                 [[277.3,-1.65],0.20],
                                 [[277.52,-2.03,0.24,0.07,120]],
                                 [[277.3,-2.07],0.13],
                                 [[278.03,-2.32,0.6,0.3,120]],
                                 [[276.5,-4.2],[276.5,-3.1],[277.8,-3.1],[277.8,-4.2],[276.5,-4.2]],
                                 ],
                     
                    },
                    
                    
                    {'name':'Aquila',      'savename':'Serpens NE',
                    'distance_error':[463,53],
                    'innerarea':[
                                 [[279,-0.5],[279,0.8],[280.2,0.8],[280.2,-0.5],[279,-0.5]],
                                 ],
                    
                    },
                    
                   {'name':'Orion A',           'savename':'Orion A',     
                    
                    'innerarea':[[[88,-11],[82,-11],[82,-3],[88,-3],[88,-11]],
                                 [[83.25,-5],[83.25,-3.75],[84.25,-3.75],[84.25,-5],[83.25,-5]],
                                 [[83.25,-5.33],[83.25,-5],[84.25,-5],[84.25,-5.33],[83.25,-5.33]],
                                 
                                 [[83.25,-5.33],[83.25,-6],[84.25,-6],[84.25,-5.33],[83.25,-5.33]],
                                 [[83.75,-7],[83.25,-6],[84.25,-6],[84.80,-7],[83.75,-7]],
                                   
                                 [[84.0,-7],[84.0,-7.5],[85.0,-8.5],[85.75,-7.9],[84.80,-7],[84.0,-7]],
                                 [[85.0,-10.2],[85.0,-8.5],[85.75,-7.9],[86.2,-8.5],[86,-10.2],[85.0,-10.2]],
                             
                                 ],
                    
                    },
                   {'name':'Auriga/CMC',      'savename':'AurigaCMC',
                    'distance_error':[496,59],
                    'innerarea':[[[68,34],[60,34],[60,42],[68,42],[68,34]],
                                 [[66.5,34.5],[66.5,37],[68,37],[68,34.5],[66.5,34.5]]],
                    
                   },
                      ]


FLUX_ALL_X = []
FLUX_ALL_Y = []
ALPHA_PC_F_select = []

alpha_all_PCs = []
for num in range(len(Protocluster)):
    coor_spitzer_distance,coor_Gaia,dis_distance,dis_low_distance,dis_up_distance,id_source_dis = get_distance(Protocluster,num)
    coor_Gaia_spitzer_pp,pm_ra_pp,pm_ra_error_pp,pm_dec_pp,pm_dec_error_pp,parallax_pp,parallax_error_pp,id_source_pp_pp = get_pm_parallax(Protocluster,num)
    coor_Gaia_spitzer,pm_ra,pm_ra_error,pm_dec,pm_dec_error,parallax,parallax_error,id_source_pp = [],[],[],[],[],[],[],[]
    dis,dis_low,dis_up = [],[],[]
    
    for i in range(len(coor_Gaia_spitzer_pp)):
        for j in range(len(coor_spitzer_distance)):
            if id_source_pp_pp[i] == id_source_dis[j]:
                coor_Gaia_spitzer.append(coor_Gaia_spitzer_pp[i])
                pm_ra.append(pm_ra_pp[i])
                pm_ra_error.append(pm_ra_error_pp[i])
                pm_dec.append(pm_dec_pp[i])
                pm_dec_error.append(pm_dec_error_pp[i])
                parallax.append(parallax_pp[i])
                parallax_error.append(parallax_error_pp[i])
                id_source_pp.append(id_source_pp_pp[i])
                
                dis.append(dis_distance[j])
                dis_low.append(dis_low_distance[j])
                dis_up.append(dis_up_distance[j])
                break

    coor_Gaia_spitzer_select,dis_select,dis_low_select,dis_up_select = [],[],[],[]
    pm_ra_select, pm_ra_error_select, pm_dec_select, pm_dec_error_select = [],[],[],[]
    
    coor_Gaia_spitzer_abandon,dis_abandon,dis_low_abandon,dis_up_abandon = [],[],[],[]
    pm_ra_abandon, pm_ra_error_abandon, pm_dec_abandon, pm_dec_error_abandon = [],[],[],[]
    
    id_source_pp_select,id_source_pp_abandon = [],[]
    for i in range(len(dis)):
        if dis_up[i] < Protocluster[num]['distance_select'][0] or \
            dis_low[i] > Protocluster[num]['distance_select'][1]: 
            coor_Gaia_spitzer_abandon.append(coor_Gaia_spitzer[i])
            dis_abandon.append(dis[i])
            dis_low_abandon.append(dis_low[i])
            dis_up_abandon.append(dis_up[i])
            
            pm_ra_abandon.append(pm_ra[i])
            pm_ra_error_abandon.append(pm_ra_error[i])
            pm_dec_abandon.append(pm_dec[i])
            pm_dec_error_abandon.append(pm_dec_error[i])
            id_source_pp_abandon.append(id_source_pp[i])
        else:
            coor_Gaia_spitzer_select.append(coor_Gaia_spitzer[i])
            dis_select.append(dis[i])
            dis_low_select.append(dis_low[i])
            dis_up_select.append(dis_up[i])
            
            pm_ra_select.append(pm_ra[i])
            pm_ra_error_select.append(pm_ra_error[i])
            pm_dec_select.append(pm_dec[i])
            pm_dec_error_select.append(pm_dec_error[i])
            id_source_pp_select.append(id_source_pp[i])
    
    coor_Gaia_spitzer_abandon_foreground = []
    coor_Gaia_spitzer_abandon_background = []
    for i in range(len(dis)):
        if dis_up[i] < Protocluster[num]['distance_select'][0]: 
            coor_Gaia_spitzer_abandon_foreground.append(coor_Gaia_spitzer[i])
        if dis_low[i] > Protocluster[num]['distance_select'][1]: 
            coor_Gaia_spitzer_abandon_background.append(coor_Gaia_spitzer[i])
    
    coor,alpha_PC = Get_coor_alpha_with_24(Protocluster,num,
                                           "Taurus_known.csv",
                                           "Taurus_new.csv",
                                           "Orion.csv",
                                           "18 protoclusters.csv")
    coor_consider,alpha_consider_with_24 = [],[]
    for i in range(len(coor)):
        if coor[i] not in coor_Gaia_spitzer_abandon:
            coor_consider.append(coor[i])
            alpha_consider_with_24.append(alpha_PC[i])
    alpha_all_PCs.append(alpha_consider_with_24)
    
    coor_F,alpha_PC_F,FLUX_x,FLUX_y = Get_coor_alpha_with_24_flux(Protocluster,num,
                                                                   "Taurus_known.csv",
                                                                   "Taurus_new.csv",
                                                                   "Orion.csv",
                                                                   "18 protoclusters.csv")
    flux_x,flux_y = [],[]
    alpha_PC_F_select = []
    for i in range(len(coor_F)):
        if coor_F[i] not in coor_Gaia_spitzer_abandon:
            flux_x.append(FLUX_x[i])
            flux_y.append(FLUX_y[i])
            alpha_PC_F_select.append(alpha_PC_F[i])
    FLUX_ALL_X.append(flux_x)
    FLUX_ALL_Y.append(flux_y)
    ALPHA_PC_F_select.append(alpha_PC_F_select)

alpha_all_in_one = []
for i in alpha_all_PCs:
    print(len(i))
    for j in i:
        alpha_all_in_one.append(j)
        

X = np.linspace(min(alpha_all_in_one),max(alpha_all_in_one),500)
counts_norm_all,bin_edges = np.histogram(alpha_all_in_one,bins='fd',density=True)
bin_all = bin_edges
bin_width = bin_edges[1]-bin_edges[0]
X = (bin_edges[1:] + bin_edges[:-1])/2
X_all = list(X)

len(alpha_all_in_one)




###revise alpha distribution of Liu et al 2023
alpha_PCs_consider_with_24 = []
PC_consider_with_24 = []

for num in range(len(Protocluster)):
    coor_spitzer_distance,coor_Gaia,dis_distance,dis_low_distance,dis_up_distance,id_source_dis = get_distance(Protocluster,num)
    coor_Gaia_spitzer_pp,pm_ra_pp,pm_ra_error_pp,pm_dec_pp,pm_dec_error_pp,parallax_pp,parallax_error_pp,id_source_pp_pp = get_pm_parallax(Protocluster,num)
    coor_Gaia_spitzer,pm_ra,pm_ra_error,pm_dec,pm_dec_error,parallax,parallax_error,id_source_pp = [],[],[],[],[],[],[],[]
    dis,dis_low,dis_up = [],[],[]
    for i in range(len(coor_Gaia_spitzer_pp)):
        for j in range(len(coor_spitzer_distance)):
            if id_source_pp_pp[i] == id_source_dis[j]:
                coor_Gaia_spitzer.append(coor_Gaia_spitzer_pp[i])
                pm_ra.append(pm_ra_pp[i])
                pm_ra_error.append(pm_ra_error_pp[i])
                pm_dec.append(pm_dec_pp[i])
                pm_dec_error.append(pm_dec_error_pp[i])
                parallax.append(parallax_pp[i])
                parallax_error.append(parallax_error_pp[i])
                id_source_pp.append(id_source_pp_pp[i])
                
                dis.append(dis_distance[j])
                dis_low.append(dis_low_distance[j])
                dis_up.append(dis_up_distance[j])
                break
            
    coor_Gaia_spitzer_select,dis_select,dis_low_select,dis_up_select = [],[],[],[]
    pm_ra_select, pm_ra_error_select, pm_dec_select, pm_dec_error_select = [],[],[],[]
    
    coor_Gaia_spitzer_abandon,dis_abandon,dis_low_abandon,dis_up_abandon = [],[],[],[]
    pm_ra_abandon, pm_ra_error_abandon, pm_dec_abandon, pm_dec_error_abandon = [],[],[],[]
    
    id_source_pp_select,id_source_pp_abandon = [],[]
    for i in range(len(dis)):
        if dis_up[i] < Protocluster[num]['distance_select'][0] or \
            dis_low[i] > Protocluster[num]['distance_select'][1]: 
            coor_Gaia_spitzer_abandon.append(coor_Gaia_spitzer[i])
            dis_abandon.append(dis[i])
            dis_low_abandon.append(dis_low[i])
            dis_up_abandon.append(dis_up[i])
            
            pm_ra_abandon.append(pm_ra[i])
            pm_ra_error_abandon.append(pm_ra_error[i])
            pm_dec_abandon.append(pm_dec[i])
            pm_dec_error_abandon.append(pm_dec_error[i])
            id_source_pp_abandon.append(id_source_pp[i])
        else:
            coor_Gaia_spitzer_select.append(coor_Gaia_spitzer[i])
            dis_select.append(dis[i])
            dis_low_select.append(dis_low[i])
            dis_up_select.append(dis_up[i])
            
            pm_ra_select.append(pm_ra[i])
            pm_ra_error_select.append(pm_ra_error[i])
            pm_dec_select.append(pm_dec[i])
            pm_dec_error_select.append(pm_dec_error[i])
            id_source_pp_select.append(id_source_pp[i])
    
    coor_Gaia_spitzer_abandon_foreground = []
    coor_Gaia_spitzer_abandon_background = []
    for i in range(len(dis)):
        if dis_up[i] < Protocluster[num]['distance_select'][0]: 
            coor_Gaia_spitzer_abandon_foreground.append(coor_Gaia_spitzer[i])
        if dis_low[i] > Protocluster[num]['distance_select'][1]: 
            coor_Gaia_spitzer_abandon_background.append(coor_Gaia_spitzer[i])
    coor,alpha_PC = Get_coor_alpha_with_24(Protocluster,num,
                                           "Taurus_known.csv",
                                           "Taurus_new.csv",
                                           "Orion.csv",
                                           "18 protoclusters.csv")
    coor_consider,alpha_consider_with_24 = [],[]
    for i in range(len(coor)):
        if coor[i] not in coor_Gaia_spitzer_abandon:
            coor_consider.append(coor[i])
            alpha_consider_with_24.append(alpha_PC[i])
            
    if len(alpha_consider_with_24) > 50:
        alpha_PCs_consider_with_24.append(alpha_consider_with_24)
        PC_consider_with_24.append(Protocluster[num]['savename'])

counts_norm_ALL = []
for i in range(len(alpha_PCs_consider_with_24)):
    counts_norm,bin_edges = np.histogram(alpha_PCs_consider_with_24[i],bins=bin_all,density=True)
    counts_norm_ALL.append(counts_norm)
# 转置二维数组
counts_norm_ALL_convert = []
for i in range(len(counts_norm_ALL[0])):  # 行数
    t = []
    for j in range(len(counts_norm_ALL)):
        t.append(counts_norm_ALL[j][i])
    counts_norm_ALL_convert.append(t)
counts_norm_average_I = [np.mean(i) for i in counts_norm_ALL_convert]
counts_norm_std = [np.std(i) for i in counts_norm_ALL_convert]

tck_PaperI = splrep(X_all,counts_norm_average_I,k=3,s=0.08)
Y_PaperI = splev(X,tck_PaperI)

All_time_I = splint(max(alpha_all_in_one),min(alpha_all_in_one),tck_PaperI)*2/splint(-1.6,-0.3,tck_PaperI)
time_I = []
for i in X:
    time_I.append(splint(5,i,tck_PaperI)*All_time_I/splint(max(alpha_all_in_one),min(alpha_all_in_one),tck_PaperI))





###new 36 sub-regions
label_part = ['A','B','C','D','E','F','G','H','I','G','K','L','M','N','O','P','Q','R','S','T','All']
alpha_all = []
Name_small_part = []
for i_care in range(len(Protocluster_cared)):#
    print(Protocluster_cared[i_care]['savename'])
    for num in range(len(Protocluster)):
        if Protocluster[num]['savename'] == Protocluster_cared[i_care]['savename']:
            
            coor_initial,alpha_initial = Get_coor_alpha_with_24_no_revise(Protocluster,num,
                                                               "Taurus_known.csv",
                                                               "Taurus_new.csv",
                                                               "Orion.csv",
                                                               "18 protoclusters.csv")
            
            coor_Gaia_distance,coor_Gaia,dis_distance,dis_low_distance,dis_up_distance,id_source_dis = get_distance(Protocluster,num)
            coor_Gaia_spitzer_pp,pm_ra_pp,pm_ra_error_pp,pm_dec_pp,pm_dec_error_pp,parallax_pp,parallax_error_pp,id_source_pp_pp = get_pm_parallax(Protocluster,num)
            coor_Gaia_spitzer,pm_ra,pm_ra_error,pm_dec,pm_dec_error,parallax,parallax_error,id_source_pp = [],[],[],[],[],[],[],[]
            dis,dis_low,dis_up = [],[],[]
            for i in range(len(coor_Gaia_spitzer_pp)):
                for j in range(len(coor_Gaia_distance)):
                    if id_source_pp_pp[i] == id_source_dis[j]:
                        coor_Gaia_spitzer.append(coor_Gaia_spitzer_pp[i])
                        pm_ra.append(pm_ra_pp[i])
                        pm_ra_error.append(pm_ra_error_pp[i])
                        pm_dec.append(pm_dec_pp[i])
                        pm_dec_error.append(pm_dec_error_pp[i])
                        parallax.append(parallax_pp[i])
                        parallax_error.append(parallax_error_pp[i])
                        id_source_pp.append(id_source_pp_pp[i])
                        
                        dis.append(dis_distance[j])
                        dis_low.append(dis_low_distance[j])
                        dis_up.append(dis_up_distance[j])
                        break
            
            #print(len(coor_Gaia_spitzer),len(pm_ra),len(pm_ra_error),len(pm_dec),len(pm_dec_error),len(parallax),len(parallax_error),len(id_source_pp),len(dis),len(dis_low),len(dis_up))
            coor_Gaia_spitzer_select,dis_select,dis_low_select,dis_up_select = [],[],[],[]
            pm_ra_select, pm_ra_error_select, pm_dec_select, pm_dec_error_select = [],[],[],[]
            
            coor_Gaia_spitzer_abandon,dis_abandon,dis_low_abandon,dis_up_abandon = [],[],[],[]
            pm_ra_abandon, pm_ra_error_abandon, pm_dec_abandon, pm_dec_error_abandon = [],[],[],[]
            
            id_source_pp_select,id_source_pp_abandon = [],[]
            for i in range(len(dis)):
                if dis_up[i] < Protocluster[num]['distance_select'][0] or \
                    dis_low[i] > Protocluster[num]['distance_select'][1]: 
                    coor_Gaia_spitzer_abandon.append(coor_Gaia_spitzer[i])
                    dis_abandon.append(dis[i])
                    dis_low_abandon.append(dis_low[i])
                    dis_up_abandon.append(dis_up[i])
                    
                    pm_ra_abandon.append(pm_ra[i])
                    pm_ra_error_abandon.append(pm_ra_error[i])
                    pm_dec_abandon.append(pm_dec[i])
                    pm_dec_error_abandon.append(pm_dec_error[i])
                    id_source_pp_abandon.append(id_source_pp[i])
                else:
                    coor_Gaia_spitzer_select.append(coor_Gaia_spitzer[i])
                    dis_select.append(dis[i])
                    dis_low_select.append(dis_low[i])
                    dis_up_select.append(dis_up[i])
                    
                    pm_ra_select.append(pm_ra[i])
                    pm_ra_error_select.append(pm_ra_error[i])
                    pm_dec_select.append(pm_dec[i])
                    pm_dec_error_select.append(pm_dec_error[i])
                    id_source_pp_select.append(id_source_pp[i])
            
            coor,alpha = [],[]
            for coor_initial_i in range(len(coor_initial)):
                if coor_initial[coor_initial_i] not in coor_Gaia_spitzer_abandon:
                    coor.append(coor_initial[coor_initial_i])
                    alpha.append(alpha_initial[coor_initial_i])
            
            for m in range(len(Protocluster_cared[i_care]['innerarea'])):
                coor_in = []
                alpha_in = []
                
                if len(Protocluster_cared[i_care]['innerarea'][m]) == 2:
                    circle = Protocluster_cared[i_care]['innerarea'][m]
                    for s in range(len(coor)):
                        if is_in_circle(coor[s],circle) == True:
                            if coor[s] not in coor_in:
                                coor_in.append(coor[s])
                                alpha_in.append(alpha[s])
                if len(Protocluster_cared[i_care]['innerarea'][m]) == 1:
                    ellipse = Protocluster_cared[i_care]['innerarea'][m][0]
                    for s in range(len(coor)):
                        if is_in_ellipse(coor[s],ellipse) == True:
                            if coor[s] not in coor_in:
                                coor_in.append(coor[s])
                                alpha_in.append(alpha[s])
                if len(Protocluster_cared[i_care]['innerarea'][m]) > 2:
                    poly = Protocluster_cared[i_care]['innerarea'][m]
                    for s in range(len(coor)):
                        if is_in_poly(coor[s],poly) == True:
                            if coor[s] not in coor_in:
                                coor_in.append(coor[s])
                                alpha_in.append(alpha[s])
                            
                if m > 0 or Protocluster_cared[i_care]['savename'] == 'Serpens NE':
                    alpha_all.append(alpha_in)
                    Name_small_part.append(Protocluster_cared[i_care]['savename']+', '+label_part[m-1])
                  # time_distributed.append(-splint(alpha_in[i],5,tck)*All_time/splint(-2,5,tck))
                    
###revised fig
counts_norm_ALL = []
for i in range(len(alpha_all)):
    counts_norm,bin_edges = np.histogram(alpha_all[i],bins=bin_all,density=True)
    counts_norm_ALL.append(counts_norm)
# 转置二维数组
counts_norm_ALL_convert = []
for i in range(len(counts_norm_ALL[0])):  # 行数
    t = []
    for j in range(len(counts_norm_ALL)):
        t.append(counts_norm_ALL[j][i])
    counts_norm_ALL_convert.append(t)
counts_norm_average = [np.mean(i) for i in counts_norm_ALL_convert]
counts_norm_quartile = [[np.percentile(i,25) for i in counts_norm_ALL_convert],
                        [np.percentile(i,75) for i in counts_norm_ALL_convert]]

tck_36 = splrep(X_all,counts_norm_average,k=3,s=0.08)
Y_36 = splev(X,tck_36)

All_time = splint(max(alpha_all_in_one),min(alpha_all_in_one),tck_36)*2/splint(-1.6,-0.3,tck_36)
time = []
alpha_test = np.linspace(min(alpha_all_in_one),max(alpha_all_in_one),500)
for i in alpha_test:
    time.append(splint(max(alpha_all_in_one),i,tck_36)*All_time/splint(max(alpha_all_in_one),min(alpha_all_in_one),tck_36))






X_Paper2023 = np.array([  -2.25,
                -1.75,
                -1.25,
                -0.75,
                -0.25,
                0.25,
                0.75,
                1.25,
                1.75,
                2.25,
                2.75,
                3.25,
                3.75,
                4.25,
                4.75])
counts_norm_average_Paper2023 = [  0,
                         0.13854081832069495,
                         0.6311583142312382,
                         0.620367815436963,
                         0.22467851853603696,
                         0.1533586943290241,
                         0.1016316606144568,
                         0.0660647230975753,
                         0.04161049874396486,
                         0.015506519925102354,
                         0.004536571730834591,
                         0.0009956080621651698,
                         0.0009324009324009324,
                         0.0,
                         0.0006178560395427864]

tck_Paper2023=splrep(X_Paper2023,counts_norm_average_Paper2023,k=2)
x_Paper2023 = np.linspace(-2,max(alpha_all_in_one),500)
y_Paper2023 = splev(x_Paper2023,tck_Paper2023)*splint(-2,max(alpha_all_in_one),tck_Paper2023)

All_time_2023 = splint(5,-2,tck_Paper2023)*2/splint(-1.6,-0.3,tck_Paper2023)
time_2023 = []
alpha_test_2023 = np.linspace(-2,max(alpha_all_in_one),500)
for i in alpha_test_2023:
    time_2023.append(splint(5,i,tck_Paper2023)*All_time_2023/splint(5,-2,tck_Paper2023))










h_fig,w_fig = 2,1
panel_size = 6
fig = plt.figure(figsize=(w_fig * panel_size ,h_fig * panel_size))

ax = fig.add_subplot(h_fig,w_fig,1)

ax.step(X_all,counts_norm_all,where='mid',c='cyan',label='Single sample',
        zorder=0)

#ax.plot(X,Y_PaperI,'r--',label='13 protoclusters',zorder=1)
ax.step(X_all,counts_norm_average_I,where='mid',c='r',label='13 protoclusters',
        zorder=1)


ax.step(X_all,counts_norm_average,where='mid',c='k',label='36 sub-regions',
        zorder=2)
ax.plot(X,Y_36,'k--',zorder=2)

ax.set_xlabel(r'$\alpha$',fontsize = 15)
ax.tick_params(axis='both',labelsize=15,colors='k')
ax.grid(False)
ax.legend(fontsize=15)
ax.minorticks_on()
#y_all = splev(x1,tck_all)
#ax.plot(x1,y_all,c='gray')


ax1 = fig.add_subplot(h_fig,w_fig,2)
ax1.plot(alpha_test,time,c='k',)
ax1.plot(X,time_I,'k--',label = 'Liu+ 2023')
ax1.plot(alpha_test_2023,time_2023,'r-',label = 'Liu et al 2023')

ax1.set_xlabel(r'$\alpha$',fontsize=15)
ax1.set_ylabel(u"-Age(Myr)",fontsize=15)
ax1.tick_params(axis='both',labelsize=15,colors='k')
ax1.minorticks_on()
ax1.vlines(x=-2,ymin=-3.1,ymax=0,colors="gray", ls="--", lw=1, )
ax1.legend(fontsize=15)
plt.savefig('alpha-age.png',bbox_inches = 'tight')
plt.show()







'''
x = np.linspace(-1.5,0.5,500000)
y = splev(x,tck_36)
for i in range(1,len(x)):
    if y[i] < y[i-1]:
        print(x[i],y[i])
        break
    
from scipy.optimize import fsolve

f_tck = lambda x: splev(x,tck_36)-0.881/2

fsolve(f_tck,[-1.1, -0.9] )

'''


        
    



























   

















def quadratic_equation(x,A,B,C):
    return A*x**2 + B*x + C

Lamda_ALL_PCs_revised,FLUX_ALL_PCs_revised = [],[]
alpha_all = []
for i in range(len(ALPHA_PC_F_select)):
    if ALPHA_PC_F_select[i] != []:
        for j in range(len(ALPHA_PC_F_select[i])):
            alpha_all.append(ALPHA_PC_F_select[i][j])
            Lamda_ALL_PCs_revised.append(FLUX_ALL_X[i][j])
            FLUX_ALL_PCs_revised.append(FLUX_ALL_Y[i][j])
            
alpha_cared_test = np.arange(-2.5,1.1,0.2)
alpha_scale = 0.1

alpha_cared = []
alpha_Bump_mean = []
Bump_fitting = [[] for i in range(len(alpha_cared_test))]
STD = []
N_cared = []
for i_alpha_range in range(len(alpha_cared_test)):
    for i in range(len(alpha_all)):
        if alpha_all[i] < alpha_cared_test[i_alpha_range]+alpha_scale and \
           alpha_all[i] > alpha_cared_test[i_alpha_range]-alpha_scale:
           X = np.array(Lamda_ALL_PCs_revised[i])
           Y = np.array(FLUX_ALL_PCs_revised[i])-FLUX_ALL_PCs_revised[i][0]
           
           popt,pcov = curve_fit(quadratic_equation, X,Y)
           A,B,C = popt
           Bump_fitting[i_alpha_range].append(A)
    if len(Bump_fitting[i_alpha_range]) >= 10:
        alpha_cared.append(alpha_cared_test[i_alpha_range])
        STD.append(np.std(Bump_fitting[i_alpha_range]))
        alpha_Bump_mean.append(np.mean(Bump_fitting[i_alpha_range]))
        N_cared.append(len(Bump_fitting[i_alpha_range]))
        #STD.append(np.percentile(Bump_fitting[i_alpha_range],75)-np.percentile(Bump_fitting[i_alpha_range],25))
        print(len(Bump_fitting[i_alpha_range]))
    
h_fig,w_fig = 1,1
fig = plt.figure(figsize=(w_fig * panel_size ,h_fig * panel_size))
ax = fig.add_subplot(h_fig,w_fig,1)

#ax.scatter(alpha_cared,STD)
#ax.plot(alpha_cared,STD)
#ax.scatter(alpha_cared,np.array(N_cared)/100)
#ax.plot(alpha_cared,np.array(N_cared)/100)
sum(N_cared)
#ax.step(alpha_cared,np.array(N_cared)/(sum(N_cared)*0.2),where='mid',label='star number')
#np.sum(np.array(N_cared)/(sum(N_cared)*0.2))

ax.step(alpha_cared,STD,where='mid',label='std of a')
ax.step(alpha_cared,alpha_Bump_mean,label='average of a')

ax.arrow(0.5,0.8,-0.23,0.3,width=0.01,head_starts_at_zero=True,head_width=0.05,color='k')
ax.arrow(0.5,0.6,-0.1,-0.3,width=0.01,head_starts_at_zero=True,head_width=0.05,color='k')
ax.text(0.5, 0.75, 'Inclination', horizontalalignment='center',
        verticalalignment='center',fontsize = 15)
ax.arrow(-1.9,1,0.2,-0.2,width=0.01,head_starts_at_zero=True,head_width=0.05,color='k')
ax.text(-1.9,1.1, 'Gaps', horizontalalignment='center',
        verticalalignment='center',fontsize = 15)

ax.set_xlabel(r'$\alpha$',fontsize = 15)
ax.set_ylabel(r'a',fontsize = 15)
#ax.set_ylabel(r'',fontsize = 15)
ax.tick_params(axis='both',labelsize=15,colors='k')
plt.legend(fontsize=15)
plt.savefig('diversity of SED_std.png',bbox_inches = 'tight',pad_inches = 0.1)
plt.show()









'''

def func_splint(x,tck,alpha_min_test,alpha_max_test):
    return_all = []
    if isinstance(x,(int,float)):
        return splint(alpha_min_test,x,tck)/splint(alpha_min_test,alpha_max_test,tck)
    else:
        for i in x:
            return_all.append(splint(alpha_min_test,i,tck)/splint(alpha_min_test,alpha_max_test,tck))
        return_all = np.array(return_all)
        return return_all


h_fig,w_fig = 6,6
panel_size = 6
fig = plt.figure(figsize=(w_fig * panel_size ,h_fig * panel_size))
alpha_all_select = [[] for i in range(len(alpha_all))]
for i in range(len(alpha_all)):
    for j in range(len(alpha_all[i])):
        if alpha_all[i][j] >= alpha_min and \
           alpha_all[i][j] <= alpha_max:
           alpha_all_select[i].append(alpha_all[i][j])
    
    ax = fig.add_subplot(h_fig,w_fig,i+1)
    ax.hist(alpha_all_select[i],bins='fd',color='red',density=True,
            histtype='step',label='N='+str(len(alpha_all_select[i])))

    ax.plot(x1,y_s,c='gray')
    ax.set_title(Name_small_part[i],fontsize=15)
    
ax.set_xlabel(r'$\alpha$',fontsize = 15)
ax.tick_params(axis='both',labelsize=15,colors='k')
ax.grid(False)
ax.legend(fontsize=15)
plt.savefig('alpha distribution.png',bbox_inches = 'tight',pad_inches = 0.1)
plt.show()

Number = []
Number_in_Class_II = []
Name_p_value = []
Name_no_random = []
alpha_min_test,alpha_max_test = -1.6,-0.3
h_fig,w_fig = 6,6
panel_size = 6
fig = plt.figure(figsize=(w_fig * panel_size ,h_fig * panel_size))
alpha_all_select = [[] for i in range(len(alpha_all))]
for i in range(len(alpha_all)):
    for j in range(len(alpha_all[i])):
        if alpha_all[i][j] >= alpha_min_test and \
           alpha_all[i][j] <= alpha_max_test:
           alpha_all_select[i].append(alpha_all[i][j])
    
    ax = fig.add_subplot(h_fig,w_fig,i+1)
    ax.hist(alpha_all_select[i],bins='fd',color='red',density=True,
            histtype='step',label='N='+str(len(alpha_all_select[i])))

    
    statistic,pvalue = scipy.stats.kstest(alpha_all_select[i],
                                          func_splint,(tck_s,alpha_min_test,alpha_max_test))
    if pvalue <= 0.05:
        Name_no_random.append(Name_small_part[i])
    Name_p_value.append(pvalue)
    
    Number.append(len(alpha_all[i]))
    Number_in_Class_II.append(len(alpha_all_select[i]))
    
    ax.set_title(Name_small_part[i]+', Pvalue = '+format(pvalue,'.3g'),fontsize=15)
    
    x1 = np.linspace(alpha_min_test,alpha_max_test,500)
    y_s = splev(x1,tck_s)/splint(alpha_min_test,alpha_max_test,tck_s)
    ax.plot(x1,y_s,c='gray')
    
    ax.legend(fontsize=15)
ax.set_xlabel(r'$\alpha$',fontsize = 15)
ax.tick_params(axis='both',labelsize=15,colors='k')
ax.grid(False)

plt.savefig('alpha distribution_Class_II.png',bbox_inches = 'tight',pad_inches = 0.1)
plt.show()


for i in Name_p_value:
    print(format(i,'3f'))
'''
