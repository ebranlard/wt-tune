import unittest
import sys
import os
import numpy as np
MyDir=os.path.dirname(__file__)
try:
    import pybra.galib as galib
except:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    import pybra.galib as galib


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestGA(unittest.TestCase):

    def test_Map(self):
        GEN_MAP=galib.GeneMap(nBases=1, name='WS', protein_ranges=[[5,15]], protein_neutr=[10])
        self.assertEqual(GEN_MAP.encode([15] ),[1.0])
        self.assertEqual(GEN_MAP.encode([5]  ),[0.0])
        self.assertEqual(GEN_MAP.decode([0.0]),[5.0])
        self.assertEqual(GEN_MAP.decode([0.5]),[10.0])
        self.assertEqual(GEN_MAP.decode([1.0]),[15.0])

        self.assertEqual(GEN_MAP.neutralProtein(),[10.0])

        CH_MAP=galib.ChromosomeMap()
        CH_MAP.add(GEN_MAP)
        CH_MAP.add(GEN_MAP)
        self.assertEqual(CH_MAP.decode([0.0,0.0]),[5.0,5.0])
        self.assertEqual(CH_MAP.encode([10.0,10.0]),[0.5,0.5])
        self.assertEqual(CH_MAP.neutralProtein(),[10,10])
        neutr=CH_MAP.neutralProtein()
        self.assertEqual(neutr,CH_MAP.decode(CH_MAP.encode(neutr)))
        
 
if __name__ == '__main__':
    unittest.main()
