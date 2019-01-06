import unittest
from unittest import skip
import time

import runfailuretest
import failurecriteria
import clt
import numpy


def assemble_valid_laminate():
    # Valid (but fictitious) material properties
    mat1 = {
     "E1"  : 10e9,       
     "E2"  : 1e9,       
     "n12" : 0.1,         
     "G12" : 1e9,       
     "Xt" : 10e6,         
     "Xc" : 10e6,
     "Yt" : 10e6,
     "Yc" : 10e6,
     "S12" : 10e6,
     "S32" : 10e6,
     "a1" : 1e-6,
     "a2" : 1e-6,
     "b1" : 0.01,
     "b2" : 0.35
     }

    # Assembles material list
    mat_list = [mat1, []]

    # Initializes dictionary of the laminate layup configurations
    # thk is thickness; ang is angle; mat_id is material id
    lam = {"thk": [], "ang": [], "mat_id" : []}

    # Adds 10 layers of 0.1mm thickness and material id 0 
    for i in range(10):
        lam["thk"].append(1e-4)
        lam["mat_id"].append(0)
    lam["ang"].extend((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

    return (lam, mat_list)

class CLTCompleteTest(unittest.TestCase):
    def test_basic_valid_laminate(self):
        initial_time = time.perf_counter()
        (lam, mat_list) = assemble_valid_laminate()
        F = [1e6, 0, 0, 0, 0, 0]
        results = clt.calc_stressCLT(mat_list, lam, F)
        finish_time = time.perf_counter()
        exec_time = finish_time - initial_time
        print("Execution time: {} seconds.".format(exec_time))

class CLTFunctionsTest(unittest.TestCase):

    def test_invalid_layup_to_assemble_Z_raises_error(self):
        """assemble_Z must not accept invalid argument."""
        lam_list = (None, '', [0, 0, 0])
        for lam in lam_list:
            self.assertRaises(clt.LaminateLayupError, clt.assemble_Z, lam)

    def test_valid_layup_to_assemble_Z_returns_valid_numpy_array(self):
        """Z_vector returned by the function must always be a np.ndarray"""
        (lam, _) = assemble_valid_laminate()
        Z_vector = clt.assemble_Z(lam)
        self.assertTrue(isinstance(Z_vector, numpy.ndarray))

    def test_known_layup_returns_expected_Z_vector(self):
        """Sends a lam with known correct Z_vector and compares to returned"""
        (lam, _) = assemble_valid_laminate()
        Z_vector = clt.assemble_Z(lam)
        Z_vector = [round(z, 6) for z in Z_vector]
        expected_Z_vector = [-5e-04, -4e-04, -3e-04, -2e-04,
                            -1e-04, 0,  1e-04,  2e-04,
                             3e-04,  4e-04,  5e-04]
        self.assertEqual(Z_vector, expected_Z_vector)

    def test_assemble_matrixT_invalid_input_returns_error(self):
        """Sends invalid inputs to the function, expects errors"""
        ang_list = (None, 361, -361, -360.1)
        for ang in ang_list:
            self.assertRaises(clt.LaminateLayupError, clt.assemble_matrixT, ang)

    def test_assemble_matrixT_valid_input_returns_numpy_array(self):
        """Sends valid input and expects a numpy array in return"""
        ang_list = (0, 90, 360, -180.9, -359.1)
        for ang in ang_list:
            T = clt.assemble_matrixT(ang)
            self.assertTrue(isinstance(T, numpy.ndarray))

    def test_assemble_matrixT_returns_expected_output_with_known_input(self):
        """Sends input that has known T matrix, expects the same output"""
        ang_list = [0, 90, -45, +45, 90, 0]
        for ang in ang_list:
            T = clt.assemble_matrixT(ang)
            T = numpy.round_(T, 6)
            if ang == 0:
                expected_T = numpy.array([[ 1,  0,  0],
                                          [ 0,  1,  0],
                                          [ 0,  0,  1]])
            elif ang == 90:
                expected_T = numpy.array([[ 0,  1,  0],
                                          [ 1,  0,  0],
                                          [ 0,  0, -1]])
            elif ang == -45:
                expected_T = numpy.array([[ 0.5,  0.5, -0.5],
                                          [ 0.5,  0.5, -0.5],
                                          [  1,    -1,    0]])
            elif ang == +45:
                expected_T = numpy.array([[ 0.5,  0.5,  0.5],
                                          [ 0.5,  0.5, -0.5],
                                          [  -1,    1,    0]])
            expected_Ti_list = [int(Ti) for Ti in numpy.nditer(expected_T)]
            returned_Ti_list = [int(Ti) for Ti in numpy.nditer(T)]
            self.assertTrue(expected_Ti_list == returned_Ti_list)

    def test_assemble_matrixQ_invalid_input_returns_error(self):
        """Sends an invalid input to the function and expects the proper error"""
        mat_prop_list = (None, [1, 1, 1], 1, 0.1)
        for mat_prop in mat_prop_list:
            self.assertRaises(
                TypeError, clt.assemble_matrixQ, mat_prop
                )

    def test_assemble_matrixQ_known_input_returns_correct_result(self):
        """Sends known input to the function and expects a known output,
        allowing for 0.1% error.
        Values come from Nasa Mechanics of Laminated Composite Plates
        Page 31 - Example 1"""
        mat1 = {
             "E1"  : 20010000.0,       
             "E2"  : 1301000.0,       
             "n12" : 0.3,         
             "G12" : 1001000.0,       
             }
        expected_Q = numpy.array([[20130785,  392656,       0],
                                  [  392656, 1308853,       0],
                                  [       0,       0, 1001000]])
        returned_Q = clt.assemble_matrixQ(mat1)

        for rQ, eQ in zip(numpy.nditer(returned_Q), 
                          numpy.nditer(expected_Q)):
            if eQ == 0 or rQ == 0:
                continue
            else:
                error = rQ/eQ
            self.assertTrue(0.999 < error < 1.001)

    def test_assemble_ABD_invalid_input_errors(self):
        """Sends invalid inputs to the function, expects the proper error"""
        valid_mat_list = [{'X':1}, ]
        (valid_lam, _) = assemble_valid_laminate()
        valid_Z = clt.assemble_Z(valid_lam)

        self.assertRaises(clt.LaminateLayupError,
                        clt.assemble_ABD,
                        valid_mat_list,
                        valid_lam,
                        None # Tests for invalid Z_vector
                        ) 

        self.assertRaises(clt.LaminateLayupError,
                        clt.assemble_ABD,
                        valid_mat_list,
                        None, # Tests for invalid lam
                        valid_Z 
                        ) 

        self.assertRaises(clt.LaminateLayupError,
                        clt.assemble_ABD,
                        None, # Tests for invalid mat_list
                        valid_mat_list, 
                        valid_Z 
                        ) 

    def test_assemble_ABD_known_input_returns_correct_result(self):
        """Sends known input to the function and expects a known output,
        allowing for 0.1% error.
        Values come from Nasa Mechanics of Laminated Composite Plates
        Page 38 - Example 2"""
        mat1 = {
             "E1"  : 20010000.0,       
             "E2"  : 1301000.0,       
             "n12" : 0.3,         
             "G12" : 1001000.0,       
             }
        mat_list = (mat1, )
        # Initializes dictionary of the laminate layup configurations
        # thk is thickness; ang is angle; mat_id is material id
        lam = {"thk": [], "ang": [], "mat_id" : []}

        # Ply 1
        lam['thk'].append(0.005)
        lam['ang'].append(45)
        lam['mat_id'].append(0)

        # Ply 2
        lam['thk'].append(0.005)
        lam['ang'].append(0)
        lam['mat_id'].append(0)

        Z = clt.assemble_Z(lam)

        expected_A = numpy.array([[133420,  24735,  23523 ],
                                  [ 24735,  39325,  23523 ],
                                  [ 23523,  23523,  30819]])
        expected_B = numpy.array([[ 169.6,  -52.02, -58.80],
                                  [-52.02,  -65.59, -58.80],
                                  [-58.80,  -58.80, -52.02]])
        expected_D = numpy.array([[1.1118, 0.2061, 0.1960],
                                  [0.2061, 0.3277, 0.1960],
                                  [0.1960, 0.1960, 0.2568]])
        expected_ABD = numpy.zeros((6,6))
        expected_ABD[:3, :3] = expected_A
        expected_ABD[:3, 3:6] = expected_ABD[3:6, :3] = expected_B
        expected_ABD[3:6, 3:6] = expected_D

        returned_ABD = clt.assemble_ABD(mat_list, lam, Z)

        for rABD, eABD in zip(numpy.nditer(returned_ABD), 
                              numpy.nditer(expected_ABD)):
            if rABD == 0 or eABD == 0:
                continue
            else:
                error = rABD/eABD
            self.assertTrue(0.999 < error < 1.001)

    def test_calc_thermal_forces_invalid_input_errors(self):
        """Sends invalid inputs to the function, expects the proper error"""
        valid_mat_list = [{'X':1}, ]
        (valid_lam, _) = assemble_valid_laminate()
        valid_Z = clt.assemble_Z(valid_lam)

        self.assertRaises(clt.LaminateLayupError,
                        clt.calc_thermal_forces,
                        valid_mat_list,
                        valid_lam,
                        None # Tests for invalid Z_vector
                        ) 

        self.assertRaises(clt.LaminateLayupError,
                        clt.calc_thermal_forces,
                        valid_mat_list,
                        None, # Tests for invalid lam
                        valid_Z 
                        ) 

        self.assertRaises(clt.LaminateLayupError,
                        clt.calc_thermal_forces,
                        None, # Tests for invalid mat_list
                        valid_mat_list, 
                        valid_Z 
                        ) 

    def test_calc_thermal_forces_known_input_returns_expected_result(self):
        """Sends known input to the function and expects a known output,
        allowing for 0.1% error.
        Values come from Nasa Mechanics of Laminated Composite Plates
        Page 51 - Example 4"""
        mat1 = {
             "E1"  : 20010000.0,       
             "E2"  : 1301000.0,       
             "n12" : 0.3,         
             "G12" : 1001000.0,   
             "a1"  : -0.04e-6,
             "a2"  : 18e-6,   
             }
        mat_list = (mat1, )
        # Initializes dictionary of the laminate layup configurations
        # thk is thickness; ang is angle; mat_id is material id
        lam = {"thk": [], "ang": [], "mat_id" : []}

        # Ply 1
        lam['thk'].append(0.005)
        lam['ang'].append(0)
        lam['mat_id'].append(0)

        # Ply 2
        lam['thk'].append(0.005)
        lam['ang'].append(45)
        lam['mat_id'].append(0)

        # Ply 3
        lam['thk'].append(0.005)
        lam['ang'].append(45)
        lam['mat_id'].append(0)

        # Ply 4
        lam['thk'].append(0.005)
        lam['ang'].append(0)
        lam['mat_id'].append(0)

        Z = clt.assemble_Z(lam)

        returned_Nt = clt.calc_thermal_forces(mat_list, lam, Z, dT=-155.6)
        expected_Nt = numpy.array([-33.57, -60.42, 12.83])

        for rNt, eNt in zip(numpy.nditer(returned_Nt), 
                              numpy.nditer(expected_Nt)):
            if rNt == 0 or eNt == 0:
                continue
            else:
                error = rNt/eNt
            self.assertTrue(0.999 < error < 1.001)

    def test_calc_moisture_forces_invalid_input_errors(self):
        """Sends invalid inputs to the function, expects the proper error"""
        (valid_lam, valid_mat_list) = assemble_valid_laminate()
        valid_Z = clt.assemble_Z(valid_lam)

        self.assertRaises(clt.LaminateLayupError,
                        clt.calc_moisture_forces,
                        valid_mat_list,
                        valid_lam,
                        None # Tests for invalid Z_vector
                        ) 

        self.assertRaises(clt.LaminateLayupError,
                        clt.calc_moisture_forces,
                        valid_mat_list,
                        None, # Tests for invalid lam
                        valid_Z 
                        ) 

        self.assertRaises(clt.LaminateLayupError,
                        clt.calc_moisture_forces,
                        None, # Tests for invalid mat_list
                        valid_mat_list, 
                        valid_Z 
                        ) 

    def test_calc_moistures_forces_known_input_returns_expected_result(self):
        """Sends known input to the function and expects a known output,
        allowing for 0.1% error.
        Values come from Nasa Mechanics of Laminated Composite Plates
        Page 53 - Example 5"""
        mat1 = {
             "E1"  : 20010000.0,       
             "E2"  : 1301000.0,       
             "n12" : 0.3,         
             "G12" : 1001000.0,   
             "b1"  : 0.01,
             "b2"  : 0.35,   
             }
        mat_list = (mat1, )
        # Initializes dictionary of the laminate layup configurations
        # thk is thickness; ang is angle; mat_id is material id
        lam = {"thk": [], "ang": [], "mat_id" : []}

        # Ply 1
        lam['thk'].append(0.005)
        lam['ang'].append(0)
        lam['mat_id'].append(0)

        # Ply 2
        lam['thk'].append(0.005)
        lam['ang'].append(45)
        lam['mat_id'].append(0)

        # Ply 3
        lam['thk'].append(0.005)
        lam['ang'].append(45)
        lam['mat_id'].append(0)

        # Ply 4
        lam['thk'].append(0.005)
        lam['ang'].append(0)
        lam['mat_id'].append(0)

        Z = clt.assemble_Z(lam)

        returned_Nm = clt.calc_moisture_forces(mat_list, lam, Z, dM=0.007)
        expected_Nm = numpy.array([51.7, 60.4, -4.3])

        for rNm, eNm in zip(numpy.nditer(returned_Nm), 
                            numpy.nditer(expected_Nm)):
            if rNm == 0 or eNm == 0:
                continue
            else:
                error = rNm/eNm
            self.assertTrue(0.999 < error < 1.001)

class FailureCriteriaTest(unittest.TestCase):

    def test_tsaiwu_2D_invalid_input_errors(self):
        (lam, mat_list) = assemble_valid_laminate()
        stress_inf = stress_sup = [None, None, None] * 10
        self.assertRaises(TypeError, 
            failurecriteria.tsaiwu_2D,
            mat_list,
            lam,
            stress_inf,
            stress_sup)

    def test_tsaiwu_2D_known_input_returns_expected_result(self):
        (lam, mat_list) = assemble_valid_laminate()
        F = [1e5, 0, 0, 0, 0, 0]
        clt_results = clt.calc_stressCLT(mat_list, lam, F)
        sfTsaiWu = failurecriteria.tsaiwu_2D(
                                mat_list, 
                                lam, 
                                clt_results["MCS"]["stress"]["inf"],
                                clt_results["MCS"]["stress"]["sup"]
        )
        self.fail('Finish this test!')


if __name__ == '__main__':
    unittest.main()