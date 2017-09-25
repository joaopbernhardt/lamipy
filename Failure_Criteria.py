"""
Functions for determining failure state of a laminate based on different criteria.


Joao Paulo Bernhardt - September 2017
"""


def tsaiwu_2D(mat_list, lam, stress_inf, stress_sup):
# FUNCTION  | tsaiwu_2D & fs_tsaiwu_2D - Calculates safety factor based on Tsai-Wu 2D Criterion.
# INPUTs	| mat_list - list with dictionaries of material properties.
# 			| lam - laminate composition - contains dictionary of thicknesses, angles and material ids.
#			| stress_inf - vector containing stresses (sig1, sig2, tau) at the bottom of a layer.
#			| stress_sup - vector containing stresses (sig1, sig2, tau) at the top of a layer.
# OUTPUTs	| fs - vector containing 2 columns: safety factor (inferior) and safety factor (superior).

    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        sig1_sup = stress_sup[0][i]
        sig2_sup = stress_sup[1][i]
        tau_sup = stress_sup[2][i]
        fs_sup = fs_tsaiwu_2D(mat_prop, sig1_sup, sig2_sup, tau_sup)
        fs["fs_sup"].append(fs_sup)

        sig1_inf = stress_inf[0][i]
        sig2_inf = stress_inf[1][i]
        tau_inf = stress_inf[2][i]
        fs_inf = fs_tsaiwu_2D(mat_prop, sig1_inf, sig2_inf, tau_inf)
        fs["fs_inf"].append(fs_inf)

    return fs


def fs_tsaiwu_2D(mat_prop, sig1, sig2, tau):

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

    return (-b + (b**2 + 4*a)**(1/2))/(2*a)

#########################################################################

def maxstress_2D(mat_list, lam, stress_inf, stress_sup):
# FUNCTION  | maxstress_2D & fs_maxstress_2D - Calculates safety factor based on Maximum Stress 2D Criterion.
# INPUTs	| mat_list - list with dictionaries of material properties.
# 			| lam - laminate composition - contains dictionary of thicknesses, angles and material ids.
#			| stress_inf - vector containing stresses (sig1, sig2, tau) at the bottom of a layer.
#			| stress_sup - vector containing stresses (sig1, sig2, tau) at the top of a layer.
# OUTPUTs	| fs - vector containing 2 columns: safety factor (inferior) and safety factor (superior).

    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        sig1_sup = stress_sup[0][i]
        sig2_sup = stress_sup[1][i]
        tau_sup = stress_sup[2][i]
        fs_sup = fs_maxstress_2D(mat_prop, sig1_sup, sig2_sup, tau_sup)
        fs["fs_sup"].append(fs_sup)

        sig1_inf = stress_inf[0][i]
        sig2_inf = stress_inf[1][i]
        tau_inf = stress_inf[2][i]
        fs_inf = fs_maxstress_2D(mat_prop, sig1_inf, sig2_inf, tau_inf)
        fs["fs_inf"].append(fs_inf)

    return fs


def fs_maxstress_2D(mat_prop, sig1, sig2, tau):


    Xt = mat_prop["Xt"]
    Xc = mat_prop["Xc"]
    Yt = mat_prop["Yt"]
    Yc = mat_prop["Yc"]
    S21 = mat_prop["S12"]

    # Verify for sig1
    if sig1 > 0:
        f_t = (sig1/Xt)
    else:
        f_t = (sig1/-Xc)

    # Verify for sig2
    if sig2 > 0:
         f_c = (sig2/Yt)
    else:
        f_c = (sig2/-Yc)

    # Verify for shear
    f_s = abs(tau)/S21

    # Result FS (1 / maximum of the 3 above)
    return 1/max(f_t, f_c, f_s)

#########################################################################

def maxstrain_2D(mat_list, lam, strain_inf, strain_sup):
# FUNCTION  | maxstrain_2D & fs_maxstrain_2D - Calculates safety factor based on Maximum Strain 2D Criterion.
# INPUTs	| mat_list - list with dictionaries of material properties.
# 			| lam - laminate composition - contains dictionary of thicknesses, angles and material ids.
#			| strain_inf - vector containing strains (eps1, eps2, gamma) at the bottom of a layer.
#			| strain_sup - vector containing strains (eps1, eps2, gamma) at the top of a layer.
# OUTPUTs	| fs - vector containing 2 columns: safety factor (inferior) and safety factor (superior).

    # Get number of layers
    num = len(lam["ang"])
    fs = {"fs_inf" : [], "fs_sup" : []}

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]

        eps1_sup = strain_sup[0][i]
        eps2_sup = strain_sup[1][i]
        gamma_sup = strain_sup[2][i]
        fs_sup = fs_maxstrain_2D(mat_prop, eps1_sup, eps2_sup, gamma_sup)
        fs["fs_sup"].append(fs_sup)

        eps1_inf = strain_inf[0][i]
        eps2_inf = strain_inf[1][i]
        gamma_inf = strain_inf[2][i]
        fs_inf = fs_maxstrain_2D(mat_prop, eps1_inf, eps2_inf, gamma_inf)
        fs["fs_inf"].append(fs_inf)

    return fs


def fs_maxstrain_2D(mat_prop, eps1, eps2, gamma):

    strainXt = mat_prop["Xt"] / mat_prop["E1"]
    strainXc = mat_prop["Xc"] / mat_prop["E1"]
    strainYt = mat_prop["Yt"] / mat_prop["E2"]
    strainYc = mat_prop["Yc"] / mat_prop["E2"]
    strainS21 = mat_prop["S12"] / mat_prop["G12"]

    # Verify for eps1
    if eps1 > 0:
        f_t = (eps1/strainXt)
    else:
        f_t = (eps1/-strainXc)

    # Verify for eps2
    if eps2 > 0:
         f_c = (eps2/strainYt)
    else:
        f_c = (eps2/-strainYc)

    # Verify for gamma
    f_s = abs(gamma)/strainS21

    # Result FS (1 / maximum of the 3 above)
    return 1/max(f_t, f_c, f_s)

#########################################################################