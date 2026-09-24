"""
lamipy project - laminated composites calculations in Python.

failurecriteria.py - Module containing functions for calculating
					 safety factors according to different criteria.

INPUTS:
lam -- laminate layup characteristics
mat_list -- list with dictionaries of material properties
stress_inf -- stresses at bottom of the layer
stress_sup -- stresses at the top of the layer
strain_inf -- strains at the bottom of the layer
strain_sup -- strains at the top of the layer

OUTPUTS:
fs -- contains safety factors coupled with failure mode

Joao Paulo Bernhardt - September 2017
"""


def tsaiwu_2D(mat_list, lam, stress_inf, stress_sup):
    """ Calculates SF and mode according to Tsai-Wu criterion 
    for the whole laminate.
    """
    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        sig1_sup = stress_sup[0][i]
        sig2_sup = stress_sup[1][i]
        tau_sup = stress_sup[2][i]
        [fs_sup, mode_sup] = fs_tsaiwu_2D(mat_prop, sig1_sup, sig2_sup, tau_sup)
        fs["fs_sup"].append([fs_sup, mode_sup])

        sig1_inf = stress_inf[0][i]
        sig2_inf = stress_inf[1][i]
        tau_inf = stress_inf[2][i]
        [fs_inf, mode_inf] = fs_tsaiwu_2D(mat_prop, sig1_inf, sig2_inf, tau_inf)
        fs["fs_inf"].append([fs_inf, mode_inf])

    return fs


def fs_tsaiwu_2D(mat_prop, sig1, sig2, tau):
    """ Calculates SF and mode according to Tsai-Wu criterion (layer-wise). """

    Xt = mat_prop["Xt"]
    Xc = mat_prop["Xc"]
    Yt = mat_prop["Yt"]
    Yc = mat_prop["Yc"]
    S21 = mat_prop["S12"]

    f11 = 1/(Xt*Xc)
    f22 = 1/(Yt*Yc)
    f12 = -1/(2*(Xt*Xc*Yt*Yc)**(1/2))
    f66 = 1/(S21**2)
    f1 = 1/Xt - 1/Xc
    f2 = 1/Yt - 1/Yc

    a = f11*sig1**2 + f22*sig2**2 + f66*tau**2 + 2*f12*sig1*sig2
    b = f1*sig1 + f2*sig2;

    sf = (-b + (b**2 + 4*a)**(1/2))/(2*a)

    # Failure mode calculations  
    H1 = abs(f1*sig1 + f11*sig1**2)
    H2 = abs(f2*sig2 + f22*sig2**2)
    H6 = abs(f66*tau**2)
 
    if max(H1,H2,H6) == H1:
        mode = "fiber"        # fiber failure
    elif max(H1,H2,H6) == H2:
        mode = "matrix"        # matrix failure
    else:
        mode = "shear"        # shear failure
    
    # Returns SF & mode    
    return [sf, mode]

#########################################################################

def maxstress_2D(mat_list, lam, stress_inf, stress_sup):
    """ Calculates SF and mode according to Maximum Stress criterion 
    for the whole laminate.
    """

    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        sig1_sup = stress_sup[0][i]
        sig2_sup = stress_sup[1][i]
        tau_sup = stress_sup[2][i]
        [fs_sup, mode_sup] = fs_maxstress_2D(mat_prop, 
        									 sig1_sup, 
        									 sig2_sup, 
        									 tau_sup)
        fs["fs_sup"].append([fs_sup, mode_sup])

        sig1_inf = stress_inf[0][i]
        sig2_inf = stress_inf[1][i]
        tau_inf = stress_inf[2][i]
        [fs_inf, mode_inf] = fs_maxstress_2D(mat_prop,
        									 sig1_inf, 
        									 sig2_inf, 
        									 tau_inf)
        fs["fs_inf"].append([fs_inf, mode_inf])

    return fs


def fs_maxstress_2D(mat_prop, sig1, sig2, tau):
    """ Calc. SF and mode according to Max. Stress criterion (layer-wise). """

    Xt = mat_prop["Xt"]
    Xc = mat_prop["Xc"]
    Yt = mat_prop["Yt"]
    Yc = mat_prop["Yc"]
    S21 = mat_prop["S12"]

    # Verify for sig1
    if sig1 > 0:
        f_1 = (sig1/Xt)
    else:
        f_1 = (sig1/-Xc)

    # Verify for sig2
    if sig2 > 0:
         f_2 = (sig2/Yt)
    else:
        f_2 = (sig2/-Yc)

    # Verify for shear
    f_s = abs(tau)/S21

    f_max = max(f_1, f_2, f_s)

    # Find failure mode
    if f_max == f_1: 
    	mode = "fiber"
    elif f_max == f_2:
    	mode = "matrix"
    else:
    	mode = "shear"

    sf = 1/f_max

    # Result FS (1 / maximum of the 3 above)
    return [sf, mode]

#########################################################################

def maxstrain_2D(mat_list, lam, strain_inf, strain_sup):
    """ Calculates SF and mode according to Maximum Strain criterion 
    for the whole laminate.
    """

    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        eps1_sup = strain_sup[0][i]
        eps2_sup = strain_sup[1][i]
        gamma_sup = strain_sup[2][i]
        [fs_sup, mode_sup] = fs_maxstrain_2D(mat_prop, 
        									 eps1_sup, 
        									 eps2_sup, 
        									 gamma_sup)
        fs["fs_sup"].append([fs_sup, mode_sup])

        eps1_inf = strain_inf[0][i]
        eps2_inf = strain_inf[1][i]
        gamma_inf = strain_inf[2][i]
        [fs_inf, mode_inf] = fs_maxstrain_2D(mat_prop, 
        									 eps1_inf, 
        									 eps2_inf, 
        									 gamma_inf)
        fs["fs_inf"].append([fs_inf, mode_inf])

    return fs


def fs_maxstrain_2D(mat_prop, eps1, eps2, gamma):
    """ Calc. SF and mode according to Max. Strain criterion (layer-wise). """

    strainXt = mat_prop["Xt"] / mat_prop["E1"]
    strainXc = mat_prop["Xc"] / mat_prop["E1"]
    strainYt = mat_prop["Yt"] / mat_prop["E2"]
    strainYc = mat_prop["Yc"] / mat_prop["E2"]
    strainS21 = mat_prop["S12"] / mat_prop["G12"]

    # Verify for eps1
    if eps1 > 0:
        f_1 = (eps1/strainXt)
    else:
        f_1 = (eps1/-strainXc)

    # Verify for eps2
    if eps2 > 0:
         f_2 = (eps2/strainYt)
    else:
        f_2 = (eps2/-strainYc)

    # Verify for gamma
    f_s = abs(gamma)/strainS21

    f_max = max(f_1, f_2, f_s)

    # Find failure mode
    if f_max == f_1: 
    	mode = "fiber"
    elif f_max == f_2:
    	mode = "matrix"
    else:
    	mode = "shear"

    sf = 1/f_max

    # Result FS (1 / maximum of the 3 above)
    return [sf, mode]

#########################################################################

def hashin_2D(mat_list, lam, stress_inf, stress_sup):
    """ Calculates SF and mode according to Hashin criterion 
    for the whole laminate.
    """

    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        sig1_sup = stress_sup[0][i]
        sig2_sup = stress_sup[1][i]
        tau_sup = stress_sup[2][i]
        [fs_sup, mode_sup] = fs_hashin_2D(mat_prop, 
        								  sig1_sup, 
        								  sig2_sup, 
        								  tau_sup)
        fs["fs_sup"].append([fs_sup, mode_sup])

        sig1_inf = stress_inf[0][i]
        sig2_inf = stress_inf[1][i]
        tau_inf = stress_inf[2][i]
        [fs_inf, mode_inf] = fs_hashin_2D(mat_prop, 
        								  sig1_inf, 
        								  sig2_inf, 
        								  tau_inf)
        fs["fs_inf"].append([fs_inf, mode_inf])

    return fs


def fs_hashin_2D(mat_prop, sig1, sig2, tau):
    """ Calc. SF and mode according to Hashin criterion (layer-wise). """

    Xt = mat_prop["Xt"]
    Xc = mat_prop["Xc"]
    Yt = mat_prop["Yt"]
    Yc = mat_prop["Yc"]
    S21 = mat_prop["S12"]

    # Verify for sig1
    if sig1 >= 0:
        f_1 = (sig1/Xt)
    else:
        f_1 = -(sig1/Xc)

    # Verify for sig2
    if sig2 >= 0:
        f_2 = ((sig2/Yt)**2 + (tau/S21)**2)**0.5
    else:
        f_2 = ((sig2/Yc)**2 + (tau/S21)**2)**0.5

    f_max = max(f_1, f_2)

    if f_max == f_1:
    	mode = "fiber"
    else:
    	mode = "matrix"

    sf = 1/f_max

    # Result FS (1 / maximum of the 3 above)
    return [sf, mode]

#########################################################################