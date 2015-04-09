import unittest
from phylogeny import _sequential_binary_partition

class TestSBP(unittest.TestCase):
    def setUp(self):
        self.tree1 = "(x,y)z;"
        self.tree2 = "((b,c)a, d)root;"
        self.tree3 = "((a,b,(c,d)e)f,(g,h)i)root;"
        self.tree4 = "((a,b)c,(f,(g,h)c)a)root;"

    def test_sbp(self):
        t = TreeNode.read(StringIO(self.tree1))
        self.assertEquals(phylogenetic_basis(t),
                          {'z': array([ 0.80442968,  0.19557032])})

        t2 = TreeNode.read(StringIO(self.tree2))
        self.assertEquals(phylogenetic_basis(t2),
                          {'a': array([ 0.57597535,  0.14002925,  0.28399541]),
                           'root': array([ 0.43595159,  0.43595159,  0.12809681])})
        t3 = TreeNode.read(StringIO(self.tree3))

        with self.assertRaises(ValueError):
            phylogenetic_basis(t3)

        t4 = TreeNode.read(StringIO(self.tree4))
        with self.assertRaises(ValueError):
            phylogenetic_basis(t4)

if __name__=="__main__":
    unittest.main()
