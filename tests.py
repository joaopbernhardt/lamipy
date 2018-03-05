import unittest
import runfailuretest
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
     "b2" : 0.35}

    # Assembles material list
    mat = [mat1, []]

    # Initializes dictionary of the laminate layup configurations
    # thk is thickness; ang is angle; mat_id is material id
    lam = {"thk": [], "ang": [], "mat_id" : []}

    # Adds 10 layers of 0.1mm thickness and material id 0 
    for i in range(10):
        lam["thk"].append(1e-3)
        lam["mat_id"].append(0)
    lam["ang"].extend((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

    return lam

class CLTFunctionsTest(unittest.TestCase):

    def test_invalid_layup_to_assemble_Z_raises_error(self):
        """assemble_Z must not accept invalid argument."""
        lam_list = (None, '', [0, 0, 0])
        for lam in lam_list:
            self.assertRaises(clt.LaminateLayupError, clt.assemble_Z, lam)

    def test_valid_layup_to_assemble_Z_returns_valid_numpy_array(self):
        """Z_vector returned by the function must always be a np.ndarray"""
        lam = assemble_valid_laminate()
        Z_vector = clt.assemble_Z(lam)
        self.assertTrue(isinstance(Z_vector, numpy.ndarray))

    def test_known_layup_returns_expected_Z_vector(self):
        """Sends a lam with known correct Z_vector and compares to returned"""
        lam = assemble_valid_laminate()
        Z_vector = clt.assemble_Z(lam)
        Z_vector = [round(z, 6) for z in Z_vector]
        expected_Z_vector = [-5e-03, -4e-03, -3e-03, -2e-03,
                            -1e-03, 0,  1e-03,  2e-03,
                             3e-03,  4e-03,  5e-03]
        self.assertEqual(Z_vector, expected_Z_vector)

if __name__ == '__main__':
    unittest.main()