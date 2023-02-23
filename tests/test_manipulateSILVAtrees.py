import unittest
from Bio import AlignIO
from ete3 import Tree
from scripts.manipulateSILVAtrees import map_names_from_taxmap, read_taxonomy_map, connect_accessions_with_nucl, remove_nodes_from_tree

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
    
    def test_full_data(self):
        self.nameMapper = map_names_from_taxmap('alns/taxmap_slv_lsu_parc_138.1.txt')
        self.taxonomyMapper = read_taxonomy_map('alns/tax_slv_lsu_138.1.tre')
        self.tree = Tree('alns/tax_slv_lsu_138.1.tre', format=1)
        self.align = AlignIO.read('alns/PTC_2063_2447_2450_2452_2453_2500_2501_2504_fullArch.sto', "stockholm")
        self.nuclOne, self.nuclTwo = 0, 9
        self.assertEqual(len(self.tree.get_leaves()), 3630)
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        self.assertEqual(len(self.tree.get_leaves()), 3630)
        self.assertEqual(self.tree.children[0].children[3].children[2].nuclOfInterest[0], 'cgacaucu')
        self.assertEqual(self.tree.children[0].children[3].children[2].name, 'Nitrososphaeria')
        namesForTruncation = ['uncultured', 'Rice Cluster I', 'J07HB67', 'J07HR59', 'J07HX64']
        truncatedTree, removed_Nodes = remove_nodes_from_tree(self.tree, namesForTruncation)
        self.assertEqual(len(truncatedTree.get_leaves()), 3558)
        self.assertEqual(len(removed_Nodes), 72)

    def test_node_deletion(self):
        self.tree = Tree('alns/tax_slv_lsu_138.1.tre', format=1)
        self.align = AlignIO.read('alns/PTC_2063_2447_2450_2452_2453_2500_2501_2504_fullArch.sto', "stockholm")
        self.nuclOne, self.nuclTwo = 0, 9
        self.assertEqual(len(self.tree.get_leaves()), 3630)
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        self.assertEqual(len(self.tree.get_leaves()), 2)

if __name__ == '__main__':
    unittest.main()