"""
    lamipy project - laminated composites calculations in Python.

    Run_FailureTest.py - Module containing functions for testing 
    composite layups.

    version_0.1: Joao Paulo Bernhardt - September 2017
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
    lam["ang"].extend((45, -45, 45, -45, 45, -45, -45, 45, -45, 45, -45, 45))

    # Vector F with the applied generalized stress (unit N/m / N.m/m)
    F = [   2277,   # Nx
            284.6,   # Ny
            0e6,   # Nxy
            0e6,   # Mx
            0e6,   # My
            0e6]   # Mxy

    # Calculates stresses and strains based on CLT.
    # res = calculated results holder;
    # LCS = Laminate System; MCS = Material Coordinate System; 
    # inf = inferior; sup = superior.
    res = CLT.calc_stressCLT(mat, lam, F)

    sfTsaiWu    = FC.tsaiwu_2D(mat, 
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

Test_A()
