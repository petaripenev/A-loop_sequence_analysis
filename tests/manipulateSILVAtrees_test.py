import unittest
from Bio import AlignIO
from ete3 import Tree
from scripts.manipulateSILVAtrees import map_names_from_taxmap, read_taxonomy_map, connect_accessions_with_nucl

class TestManipulateSILVAtrees(unittest.TestCase):
    
    def setUp(self):
        self.nameMapper = map_names_from_taxmap('tests/test_data/taxmap_test.txt')
        self.taxonomyMapper = read_taxonomy_map('tests/test_data/test_tree.tre')
        self.tree = Tree('tests/test_data/test_tree.tre', format=1)
        self.align = AlignIO.read('tests/test_data/test_aln.sto', "stockholm")
        self.nuclOne, self.nuclTwo = 0, 1

    def test_map_names_from_taxmap(self):
        self.assertIsInstance(self.nameMapper, dict)
        self.assertEqual(len(self.nameMapper), 3)

    def test_read_taxonomy_map(self):
        self.assertIsInstance(self.taxonomyMapper, dict)
        self.assertEqual(len(self.taxonomyMapper), 3)

    def test_connect_accessions_with_nucl(self):
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        self.assertEqual(self.tree.children[0].nuclOfInterest[0], 'g')
        self.assertEqual(self.tree.children[1].nuclOfInterest[0], 'c')
        self.assertEqual(self.tree.children[2].nuclOfInterest[0], 'a')
        self.assertIsInstance(accessionToNucl, dict)
        self.assertEqual(len(accessionToNucl), 3)