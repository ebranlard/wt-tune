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

        
        GEN_MAP2=galib.GeneMap(nBases=2, name='RPM-PITCH', protein_ranges=[[15,25],[0,5]], protein_neutr=[20,0])
        self.assertEqual(GEN_MAP2.geneBounds(),[(0,1),(0,1)])
        self.assertEqual(GEN_MAP2.proteinBounds(),[(15,25),(0,5)])

        CH_MAP=galib.ChromosomeMap()
        CH_MAP.add(GEN_MAP)
        CH_MAP.add(GEN_MAP2)
        self.assertEqual(CH_MAP.chromosomeBounds(),[(0,1),(0,1),(0,1)])
        self.assertEqual(CH_MAP.proteinChainBounds(),[(5,15),(15,25),(0,5)])
 
        self.assertEqual(GEN_MAP2.decode(0.5,iBase = 0), 20)
        self.assertEqual(GEN_MAP2.decode(0.5,iBase = 1), 2.5)
        self.assertEqual(CH_MAP.decode  (0.5,iBase = 0), 10)
        self.assertEqual(CH_MAP.decode  (0.5,iBase = 1), 20)
        self.assertEqual(CH_MAP.decode  (0.5,iBase = 2), 2.5)

        #GEN_MAP=galib.GeneMap(nBases=1, name='WS-2', protein_ranges=[[0,100]], resolution=1000)
        #print(GEN_MAP.show_full([0]))
        #print(GEN_MAP.show_full([0.001]))
        #print(GEN_MAP.show_full([0.5]))
        #print(GEN_MAP.show_full([1]))



if __name__ == '__main__':
    unittest.main()
