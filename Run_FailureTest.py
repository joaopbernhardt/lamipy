"""
    lamipy project - laminated composites calculations in Python.

    Run_FailureTest.py - Module containing functions for testing 
    composite layups.

    Joao Paulo Bernhardt - September 2017
"""

import numpy as np
import CLT
import Failure_Criteria as FC

def Test_A():
# This temporary function tests the implementation.
    
    # Material 1 - Dictionary of properties
    mat1 = {
     "E1"  : 156.4e9,
     "E2"  : 7.786e9,
     "n12" : 0.354,
     "G12" : 3.762e9,
     "Xt" : 1826e6,
     "Xc" : 1134e6,
     "Yt" : 19e6,
     "Yc" : 131e6,
     "S12" : 75e6,
     "S32" : 41e6 }

    # Assembles material list
    mat = [mat1, []]

    # Initializes dictionary of the laminate layup configurations
    # thk = thickness; ang = angle; mat_id = material id
    lam = {"thk": [], "ang": [], "mat_id" : []}

    # Adds 12 layers of 0.127mm thickness and material id 0
    for i in range(12):
        lam["thk"].append(0.127e-3)
        lam["mat_id"].append(0)
    #lam["ang"].extend((45, -45, 45, -45, 45, -45, -45, 45, -45, 45, -45, 45))
    lam["ang"].extend((0, 0, 45, -45, 45, 90, 90, 45, -45, 45, 0, 0))

    # Vector F with the applied generalized stress (unit N/m / N.m/m)
    F = np.array([  1e4,   # Nx
                    0e4,   # Ny
                    0e4,   # Nxy
                    0e6,   # Mx
                    0e6,   # My
                    0e4])   # Mxy

    # Calculates stresses and strains based on CLT.
    # res = calculated results holder;
    # LCS = Laminate System; MCS = Material Coordinate System; 
    # inf = inferior; sup = superior.
    res = CLT.calc_stressCLT(mat, lam, F)

    sfTsaiWu    = FC.tsaiwu_2D(   mat, 
                                  lam, 
                                  res["MCS"]["stress"]["inf"],
                                  res["MCS"]["stress"]["sup"])
    
    sfMaxStress = FC.maxstress_2D(mat, 
                                  lam, 
                                  res["MCS"]["stress"]["inf"], 
                                  res["MCS"]["stress"]["sup"])
    
    sfMaxStrain = FC.maxstrain_2D(mat, 
                                  lam, 
                                  res["MCS"]["strain"]["inf"], 
                                  res["MCS"]["strain"]["sup"])

    sfHashin = FC.hashin_2D(      mat, 
                                  lam, 
                                  res["MCS"]["stress"]["inf"], 
                                  res["MCS"]["stress"]["sup"])

    Progressive_Failure_Test(mat, lam, F)

    pass

def Progressive_Failure_Test(mat, lam, F):

    num = len(lam["ang"])
    #failed_list = [[False, ""]] * num
    fail_status =  {"Failed?" : [False] * num, 
                    "Mode" : [""] * num,
                    "Load Factor" : [0] * num}
    failed_count = [0, 0]
    LF = 1.00       # Load factor
    LS = 1.05       # Load step multiplier

    while failed_count[0] < num:

        # Sends for CLT calculations
        res = CLT.calc_stressCLT(mat, lam, F*LF, fail_status["Mode"])

        # Sends for SF calculations
        SF_list = FC.tsaiwu_2D(   mat, 
                                  lam, 
                                  res["MCS"]["stress"]["inf"],
                                  res["MCS"]["stress"]["sup"])

        # Loop for layer failure check
        for i in range(num):
            
            # Updates sf values
            sf_inf = SF_list["fs_inf"][i][0]
            sf_sup = SF_list["fs_sup"][i][0]
            sf = min(sf_inf, sf_sup)
            
            # Gets failure mode based on min sf
            if sf == SF_list["fs_inf"][i][0]:
                mode = SF_list["fs_inf"][i][1]
            else:
                mode = SF_list["fs_sup"][i][1]

            # Checks for new failure
            if sf < 1 and not fail_status["Failed?"][i]:
                fail_status["Failed?"][i] = True
                fail_status["Mode"][i] = mode
                fail_status["Load Factor"][i] = LF
                failed_count[1] = failed_count[1] + 1
                #print("Layer "+str(i)+" has failed. Mode: " + mode)

        # Increases LF if no new failure
        if failed_count[1] == failed_count[0]:       
            LF = LF*LS        #increases Load Factor by 5%
            #print([int(load) for load in LF*F if load>0])

        failed_count[0] = failed_count[1]

    fpf = min(fail_status["Load Factor"])
    lpf = max(fail_status["Load Factor"])

    print("First Ply Failure at LF: " + str(round(fpf)))
    print("Last Ply Failure at LF: " + str(round(lpf)))
    print("LPF / FPF : " + str(round(lpf/fpf, 1)))
    
    pass

Test_A()