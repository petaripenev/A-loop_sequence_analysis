import unittest
import matplotlib.cm as cm
from Bio import AlignIO
from ete3 import Tree
from scripts.manipulateSILVAtrees import map_names_from_taxmap, read_taxonomy_map
from scripts.manipulateSILVAtrees import connect_accessions_with_nucl, remove_nodes_from_tree
from scripts.manipulateSILVAtrees import combine_tempura_with_tree, save_colored_tree_by_tempura
from scripts.manipulateSILVAtrees import generate_nucleotide_type_set, save_shape_annotations

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
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, 
            self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        self.assertEqual(self.tree.children[0].nuclOfInterest[0], 'g')
        self.assertEqual(self.tree.children[1].nuclOfInterest[0], 'c')
        self.assertEqual(self.tree.children[2].nuclOfInterest[0], 'a')
        self.assertIsInstance(accessionToNucl, dict)
        self.assertEqual(len(accessionToNucl), 3)
    
    def test_full_data(self):
        self.nameMapper = map_names_from_taxmap('alns/taxmap_slv_lsu_parc_138.1.txt')
        self.taxonomyMapper = read_taxonomy_map('alns/tax_slv_lsu_138.1.tre')
        self.tree = Tree('alns/tax_slv_lsu_138.1.tre', format=1)
        self.align = AlignIO.read('alns/PTC_2063_2447_2450_2452_2453_2500_2501_2504_fullArch.sto', 
            "stockholm")
        self.nuclOne, self.nuclTwo = 0, 9
        self.assertEqual(len(self.tree.get_leaves()), 3630)
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, 
            self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        self.assertEqual(len(self.tree.get_leaves()), 3630)
        self.assertEqual(self.tree.children[0].children[3].children[2].nuclOfInterest[0], 
            'cgacaucu')
        self.assertEqual(self.tree.children[0].children[3].children[2].name, 
            'Nitrososphaeria')
        namesForTruncation = ['uncultured', 'Rice Cluster I', 'J07HB67', 'J07HR59', 'J07HX64']
        truncatedTree, removed_Nodes = remove_nodes_from_tree(self.tree, namesForTruncation)
        self.assertEqual(len(truncatedTree.get_leaves()), 3558)
        self.assertEqual(len(removed_Nodes), 72)

    def test_node_deletion(self):
        self.tree = Tree('alns/tax_slv_lsu_138.1.tre', format=1)
        self.align = AlignIO.read('alns/PTC_2063_2447_2450_2452_2453_2500_2501_2504_fullArch.sto', 
            "stockholm")
        self.nuclOne, self.nuclTwo = 0, 9
        self.assertEqual(len(self.tree.get_leaves()), 3630)
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, 
            self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        self.assertEqual(len(self.tree.get_leaves()), 2)
    
    def test_combine_tempura_with_tree(self):
        self.nameMapper = map_names_from_taxmap('alns/taxmap_slv_lsu_parc_138.1.txt')
        self.taxonomyMapper = read_taxonomy_map('alns/tax_slv_lsu_138.1.tre')
        self.tree = Tree('alns/tax_slv_lsu_138.1.tre', format=1)
        self.align = AlignIO.read('alns/PTC_2063_2447_2450_2452_2453_2500_2501_2504_fullArch.sto', 
            "stockholm")
        self.nuclOne, self.nuclTwo = 0, 9
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, 
            self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        namesForTruncation = ['uncultured', 'Rice Cluster I', 'J07HB67', 'J07HR59', 'J07HX64']
        truncatedTree, removed_Nodes = remove_nodes_from_tree(self.tree, namesForTruncation)
        toptList, toptMatchList, tempLabeledLeaves = combine_tempura_with_tree('tests/test_data/TEMPURA_test.csv', 
            truncatedTree)
        self.assertEqual(len(toptList), 5)
        self.assertEqual(len(toptMatchList), 3)
        self.assertEqual(len(tempLabeledLeaves), 3)
        self.assertEqual(toptList[0], 100.0)
        self.assertEqual(toptMatchList[1], 105.0)
        self.assertEqual(tempLabeledLeaves[0], 'Methanopyrus')
        has_temp = [hasattr(x,'topt_ave') for x in truncatedTree.traverse() 
            if x.name == tempLabeledLeaves[0]]
        self.assertEqual(has_temp[0], True)
        mapper = cm.ScalarMappable()
        save_colored_tree_by_tempura('test', truncatedTree, toptList, mapper, 
            location='./tests/test_data/')
        with open('tests/test_data/temperatureLabelColors_test.txt', 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertEqual(line, 'DATASET_STYLE\n')
                elif i == 1:
                    self.assertEqual(line, 'SEPARATOR COMMA\n')
                elif i == 2:
                    self.assertEqual(line, 'DATASET_LABEL,Colors 30.0-106.0\n')
                elif i == 3:
                    self.assertEqual(line, 'COLOR,#ffff00\n')
                elif i == 4:
                    self.assertEqual(line, 'DATA\n')
                elif i == 53:
                    self.assertEqual(line.split(',')[0], 'Methanopyrus')
                else:
                    self.assertEqual(len(line.split(',')), 7)
            self.assertEqual(i, 3562)

    def test_save_shape_annotations(self):
        self.nameMapper = map_names_from_taxmap('alns/taxmap_slv_lsu_parc_138.1.txt')
        self.taxonomyMapper = read_taxonomy_map('alns/tax_slv_lsu_138.1.tre')
        self.tree = Tree('alns/tax_slv_lsu_138.1.tre', format=1)
        self.align = AlignIO.read('alns/PTC_2063_2447_2450_2452_2453_2500_2501_2504_fullArch.sto', 
            "stockholm")
        self.nuclOne, self.nuclTwo = 0, 2
        accessionToNucl = connect_accessions_with_nucl(self.tree, self.align, self.nuclOne, 
            self.nuclTwo, self.nameMapper, self.taxonomyMapper)
        namesForTruncation = ['uncultured', 'Rice Cluster I', 'J07HB67', 'J07HR59', 'J07HX64']
        truncatedTree, removed_Nodes = remove_nodes_from_tree(self.tree, namesForTruncation)
        
        setOfNuclTypes = generate_nucleotide_type_set(accessionToNucl)
        self.assertEqual(len(setOfNuclTypes), 14)
        shapes = ["1","2","3","4","5","6","HH","HV","EL","DI","PL","PR","PU","PD","OC","GP"]

        truncatedTree = save_shape_annotations(shapes, setOfNuclTypes, accessionToNucl, 
            'Test', truncatedTree, outpath='tests/test_data/')
        with open(file='tests/test_data/interestNuclShapes_Test.txt', mode='r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertEqual(line, 'DATASET_SYMBOL\n')
                elif i == 1:
                    self.assertEqual(line, 'SEPARATOR COMMA\n')
                elif i == 2:
                    self.assertEqual(line, 'DATASET_LABEL,SILVA tree shapes\n')
                elif i == 3:
                    self.assertEqual(line, 'COLOR,#ff0000\n')
                elif i == 4:
                    self.assertEqual(line, 'LEGEND_TITLE,Nucleotide types,\n')
                elif i == 5:
                    self.assertEqual(line, 'LEGEND_SHAPES,1,2,3,4,5,6,HH,HV,EL,DI,PL,PR,PU,PD\n')
                elif i == 6:
                    self.assertEqual(line, 'LEGEND_COLORS,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000,#000000            \n')
                elif i == 7:
                    self.assertEqual(line, 'LEGEND_LABELS,cg,ca,-g,cc,aa,c-,-u,ag,ua,u-,uu,-a,gg,XX\n')
                elif i == 8:
                    self.assertEqual(line, 'DATA\n')
                elif i == 57:
                    self.assertEqual(line.split(',')[0], 'Methanopyrus')
                    self.assertEqual(line.split(',')[1], '1')
                else:
                    self.assertEqual(len(line.split(',')), 6)


if __name__ == '__main__':
    unittest.main()