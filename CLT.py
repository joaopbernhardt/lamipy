"""
    lamipy project - laminated composites calculations in Python.

    CLT.py - Module containing functions for stress&strain calculations 
    according to the CLASSICAL LAMINATE THEORY. 

    version_0.1: Joao Paulo Bernhardt - September 2017
"""

import numpy

def calc_stressCLT(mat_list, lam, F):
# FUNCTION  | calc_stressCLT
#           | Calculate stress and strain vectors according to the  
#           | Classical Laminate Theory, at the top and bottom of each layer.
# INPUTs	| mat_list - list with dictionaries of material properties.
#           | lam - laminate composition - contains dictionary of 
#           |       thicknesses, angles and material ids.
#			|  F - 'Force' vector (generalized stress vector) 
#           |       contains applied loads on the composite.
# OUTPUTs	| dictionary with stresses (on the material coord. system) 
#           | and strains (material & laminate coord. system.)

    # Get number of layers
    num = len(lam["ang"])
    total_thk = numpy.sum(lam["thk"])

    # Z vector contains z coordinates.
    Z = numpy.zeros(num + 1)
    Z[0] = -total_thk/2

    for i in range(num):
        Z[i + 1] = Z[i] + lam["thk"][i]

    # Calculates stiffness matrix based on laminate properties.
    ABD = assemble_ABD(mat_list, lam, Z)

    # Calculates strain vector by solving eq. form AX = B
    strain_vector = numpy.linalg.solve(ABD, F)

    # Initializes strain and curvature vectors & defines values
    strains = curvatures = numpy.zeros((1, 3))
    strains = strain_vector[:3]
    curvatures = strain_vector [3:6]

    # Initializes Laminate System (LS) strain vectors (inferior and superior)
    LS_strain_inf = numpy.zeros((3, num))
    LS_strain_sup = numpy.zeros((3, num))

    # Assign Laminate System strain (Epsilon = Epsilon0 + kappa*z)
    for i in range(num):
        LS_strain_inf[:3,i] = strains + curvatures*Z[i]
        LS_strain_sup[:3,i] = strains + curvatures*Z[i+1]

    # Initialize Material System stress & strain vectors (inferior and superior)
    MS_strain_inf = numpy.zeros((3, num))
    MS_strain_sup = numpy.zeros((3, num))
    MS_stress_inf = numpy.zeros((3, num))
    MS_stress_sup = numpy.zeros((3, num))

    # Calculates Material System stresses and strains
    for i in range(num):
        mat_id = lam["mat_id"][i]
        Q = assemble_matrixQ(mat_list[mat_id])
        T = assemble_matrixT(lam["ang"][i])
        MS_strain_sup[:,i] = numpy.dot(T, LS_strain_sup[:,i])
        MS_strain_inf[:,i] = numpy.dot(T, LS_strain_inf[:,i])
        MS_stress_sup[:,i] = numpy.dot(Q, MS_strain_sup[:,i])
        MS_stress_inf[:,i] = numpy.dot(Q, MS_strain_inf[:,i])

    # Outputs the stresses and strains vectors.
    # LCS = Laminate System; MCS = Material Coordinate System; 
    # inf = inferior; sup = superior.
    return {"LCS" : {"stress" : {"inf" : LS_strain_inf,
                                 "sup" : LS_strain_sup}}, 
            "MCS" : {"stress" : {"inf" : MS_stress_inf,
                                 "sup" : MS_stress_sup}, 
                     "strain" : {"inf" : MS_strain_inf,
                                 "sup" : MS_strain_sup}}}


##### BELOW: Auxiliary functions for the calculation of stresses and strains

def assemble_matrixQ (mat_prop):
# FUNCTION  | assemble_matrixQ
# INPUTs	| mat_prop - material properties arranged in a dictionary.
# OUTPUTs	| Q - reduced elastic matrix

    n21 = mat_prop["n12"]*mat_prop["E2"]/mat_prop["E1"]
    
    Q11 = mat_prop["E1"]/(1 - mat_prop["n12"]*n21)    
    # The Q11 equation below produces less precise results (~1% difference)
    #Q11 = mat_prop["E1"]**2 / (mat_prop["E1"] - mat_prop["n12"] * mat_prop["E2"])
    
    Q12 = (mat_prop["n12"]*mat_prop["E1"]*mat_prop["E2"] 
          / (mat_prop["E1"] - (mat_prop["n12"] ** 2) * mat_prop["E2"]))
    Q22 = (mat_prop["E1"]*mat_prop["E2"] / (mat_prop["E1"] 
          - (mat_prop["n12"] ** 2)*mat_prop["E2"]))
    Q66 = mat_prop["G12"]

    Q = numpy.zeros((3, 3))
    Q = [[Q11, Q12, 0],
         [Q12, Q22, 0],
         [0,   0, Q66]]
    return Q


def assemble_matrixT(angle):
# FUNCTION  | assemble_matrixT
#           | Transforms from Laminate Coord. Sys. to Material C. Sys.
# INPUTs	| angle (degrees) - layup angle of the layer
# OUTPUTs	| T - transformation matrix from

    #Transforms angle (degrees) to angle (radians)
    angle = numpy.pi*angle/180

    cos = numpy.cos(angle)
    sin = numpy.sin(angle)
    cs = cos*sin
    cc = cos**2
    ss = sin**2

    T = numpy.zeros((3, 3))
    T = [[cc,    ss,   cs   ],
         [ss,    cc,  -cs   ],
         [-2*cs, 2*cs, cc-ss]]
    return T


def assemble_ABD(mat_list, lam, Z):
# FUNCTION  | assemble_ABD
# INPUTs	| mat_list - list with dictionaries of material properties.
# 			| lam - laminate composition - contains dictionary of 
#           |       thicknesses, angles and material ids.
#			| Z - vector containing z-coordinates
# OUTPUTs	| ABD - laminate stiffness matrix

    num = len(lam["ang"])
    A = B = D = numpy.zeros((3,3))

    for i in range(num):
        mat_id = lam["mat_id"][i]
        mat_prop = mat_list[mat_id]
        ang = lam["ang"][i]        

        T = assemble_matrixT(ang)
        Q = assemble_matrixQ(mat_prop)
        Ti = numpy.transpose(T)
        
        #Qi = Ti*Q*T
        Qi = numpy.matmul(numpy.matmul(Ti, Q), T)

        A = A + Qi*(Z[i+1] - Z[i])
        B = B + (1/2)*Qi*(Z[i+1]**2 - Z[i]**2)
        D = D + (1/3)*Qi*(Z[i+1]**3 - Z[i]**3)

    ABD = numpy.zeros((6,6))
    ABD[:3, :3] = A
    ABD[:3, 3:6] = ABD[3:6, :3] = B
    ABD[3:6, 3:6] = D

    return ABD