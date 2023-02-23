#!/usr/bin/env python3
'''Given a SILVA .tre and map files from https://www.arb-silva.de/no_cache/download/archive/release_138.1/Exports/taxonomy/
as well as TEMPURA csv (from http://togodb.org/db/tempura), generates a color file for ITOL and adjusted tree with names and if needed truncations.'''
import csv, sys, argparse
from ete3 import NCBITaxa, Tree
from Bio import AlignIO
import matplotlib.cm as cm
from matplotlib.colors import Normalize, rgb2hex

from scripts.analyzeStockholm import organizeBPsBYindex

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('treeFile', help='Path to SILVA .tre file with phylogeny using SILVA identifiers.')
    parser.add_argument('mapFile', help='Path to SILVA taxmap file matching phylogeny names and sequence accessions to SILVA identifiers.')
    parser.add_argument('tempuraFile', help='Path to TEMPURA file with temperature data.')
    parser.add_argument('alignmentFile', help='Path to STOCKHOLM file with secondary structure data.')
    parser.add_argument('nucleotideIndex', nargs=2, help='Takes 2 integers, both indexing the nucleotide columns of interest in the STOCKHOLM file.')
    parser.add_argument('-t','--taxonomyName', help='Taxonomy name to filter the SILVA tree by (include it and all of its children: e.g. Bacteria)\n\
    Should be available in the map file.', default=None)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def map_names_from_taxmap(mapFile):
    nameMapper = dict()
    with open(mapFile) as mapFile:                                                                                          
        reader = csv.DictReader(mapFile, delimiter='\t')
        #next(reader)
        for entry in reader:
            if entry['taxid'] not in nameMapper.keys():
                nameMapper[entry['taxid']] = list()
            nameMapper[entry['taxid']].append((entry['organism_name'], entry['primaryAccession'], entry['path']))
    return nameMapper

def read_taxonomy_map(treeFile):
    taxonomyMapper = dict()
    with open(treeFile.replace('.tre','.map')) as mapFile:                                                                                          
        reader = csv.reader(mapFile, delimiter='\t')
        for entry in reader:
            taxonomyMapper[entry[0]] = entry[1]
    return taxonomyMapper

def connect_accessions_with_nucl(tree, align, nuclOne, nuclTwo, nameMapper, taxonomyMapper):
    accessionToNucl = dict()
    for node in tree.traverse("postorder"):
        if node.name in nameMapper.keys():
            accessionsFromMap = [x[1] for x in nameMapper[node.name]]
            accessionsFromAln = [x.name.split('.')[0] for x in align]
            for accession in list(set(accessionsFromMap).intersection(accessionsFromAln)):
                matchingIndexes = [i for i, item in enumerate(accessionsFromAln) if item == accession]
                for index in matchingIndexes:
                    if accession not in accessionToNucl.keys():
                        accessionToNucl[accession] = list()
                    node.add_features(nuclOfInterest=list())
                    node.add_features(sequenceNames=list())
                    #Be careful with this line, it is a bit of a hack to get the nucleotide sequence from the alignment
                    #it won't work for proper ranges of sequences, but it will work for the current case
                    accessionToNucl[accession].append(str(align[index].seq[nuclOne:nuclTwo]))
                    node.nuclOfInterest.append(str(align[index].seq[nuclOne:nuclTwo]))
                    node.sequenceNames.append(align[index].name)
            node.add_features(silvaID=node.name)
            #We want to use the name from .map file, not the ID from .tre file
            node.name = taxonomyMapper[node.name]
        elif node.name in taxonomyMapper.keys():
            node.name = taxonomyMapper[node.name]
        else:
            node.delete()
    return accessionToNucl

def remove_nodes_from_tree(truncatedTree, namesForTruncation):
    removed_Nodes = list()
    for node in truncatedTree.get_leaves():
        if node.name not in namesForTruncation:
            continue
        removedNode = node.detach()
        removed_Nodes.append(removedNode)
    return truncatedTree, removed_Nodes

def main(commandline_args):
    comm_args = create_and_parse_argument_options(commandline_args)
    tree = Tree(comm_args.treeFile, format=1)
    align = AlignIO.read(comm_args.alignmentFile, "stockholm")
    nuclOne, nuclTwo = int(comm_args.nucleotideIndex[0]), int(comm_args.nucleotideIndex[1])

    
    nameMapper = map_names_from_taxmap(comm_args.mapFile) 
    taxonomyMapper = read_taxonomy_map(comm_args.treeFile)
    
    accessionToNucl = connect_accessions_with_nucl(tree, align, nuclOne, nuclTwo, nameMapper, taxonomyMapper)

    if comm_args.taxonomyName:
        nodeOfInterest = tree.search_nodes(name=comm_args.taxonomyName)[0]
        truncatedTree = nodeOfInterest.detach()
    else:
        truncatedTree = tree   

    namesForTruncation = ['uncultured', 'Rice Cluster I', 'J07HB67', 'J07HR59', 'J07HX64']
    truncatedTree, removed_Nodes = remove_nodes_from_tree(truncatedTree, namesForTruncation)

    toptMatchList, toptList, tempLabeledLeaves = list(), list(), list()
    with open(comm_args.tempuraFile) as tempuraFile:
        next(tempuraFile)
        reader = csv.reader(tempuraFile, delimiter=',')
        for entry in reader:
            toptList.append(float(entry[15]))
            matchingNames = list(filter(lambda n: n.name.replace('Candidatus ','') in entry[0], truncatedTree.traverse() ))
            if len(matchingNames) > 0:
                matchingNames[0].add_features(topt_ave=float(entry[15]))
                tempLabeledLeaves.append(matchingNames[0].name)
                toptMatchList.append(float(entry[15]))

    #Hardcode the min to push 25 degrees to be green
    norm = Normalize(vmin=min(toptList), vmax=max(toptList), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)

    #with open(f'./temperatureTrees/temperatureLabelColors_{comm_args.taxonomyName}.txt', 'w') as f:
    #   f.write(f'DATASET_STYLE\r\nSEPARATOR COMMA\r\nDATASET_LABEL,Colors {min(toptList)}-{max(toptList)}\r\nCOLOR,#ffff00\r\nDATA\r\n')
    #   for node in truncatedTree.get_leaves():
    #       labelColor = '#ffffff'
    #       labelTextColor = '#000000'
    #       if 'topt_ave' in node.features:
    #           labelColor = rgb2hex(mapper.to_rgba(node.topt_ave))
    #           if node.topt_ave < 50:
    #               labelTextColor = '#ffffff'
    #       f.write(f'{node.name.replace("Candidatus ", "Cand. ")},label,node,{labelTextColor},1,normal,{labelColor}\r\n')
    
    downsetted = [list(set(x)) for x in accessionToNucl.values()]
    #setOfNuclTypes = sorted(list(set([list(x)[0] for x in downsetted if len(x) == 1])),reverse=True)
    only_single_entries = [list(x)[0] for x in downsetted if len(x) == 1]
    orderedNuclTypesByCount = sorted(only_single_entries, key=only_single_entries.count,reverse=True)
    setOfNuclTypes = list(dict.fromkeys(orderedNuclTypesByCount))
    setOfNuclTypes.append('XX')
    #Hardcoding for now; for later use a % of the unique nucleotide types
    #len([x[0] for x in downsetted if x[0] == 'UU'])/len(downsetted)
    #setOfNuclTypes = ['UU', 'CU', 'CC', 'UA', 'UC', 'XX']
    shapes = ["1","2","3","4","5","6","HH","HV","EL","DI","PL","PR","PU","PD","OC","GP"]
    nuclShapes = [shapes[x] if len(shapes) > x else "GP" for x in range(len(setOfNuclTypes))]
    nuclTypesToShapes = {setOfNuclTypes[i]: nuclShapes[i] for i in range(len(setOfNuclTypes))}
    stringOfNuclTypes = ','.join(nuclShapes)
    stringOfNuclColors = ','.join(['#000000' for x in range(1,len(setOfNuclTypes)+1)])
    multiMatches = [(k,v) for k,v in accessionToNucl.items() if len(v) > 1]
    with open(f'./temperatureTrees/bridge_interestNuclShapes_{comm_args.taxonomyName}.txt', 'w') as f:
        size = '15'
        f.write('DATASET_SYMBOL\r\nSEPARATOR COMMA\r\nDATASET_LABEL,SILVA tree shapes\r\nCOLOR,#ff0000\r\n')
        f.write(f'LEGEND_TITLE,Nucleotide types,\r\nLEGEND_SHAPES,{stringOfNuclTypes}\r\nLEGEND_COLORS,{stringOfNuclColors}\
            \r\nLEGEND_LABELS,{",".join(setOfNuclTypes)}\r\nDATA\r\n')
        for node in truncatedTree.get_leaves():
            labelColor = '#000000'
            if 'nuclOfInterest' in node.features:
                if len(set(node.nuclOfInterest)) > 1:
                    shape = nuclTypesToShapes['XX']
                elif node.nuclOfInterest[0] in setOfNuclTypes:
                    shape = nuclTypesToShapes[node.nuclOfInterest[0]]
                else:
                    shape = nuclTypesToShapes['XX']
            f.write(f'{node.name.replace("Candidatus ", "Cand. ")},{shape},{size},#000000,1,1\r\n')

    truncatedTree.write(format=1, outfile=f'./temperatureTrees/bridge_truncatedTree_{comm_args.taxonomyName}.nwk')




if __name__ == "__main__":
    main(sys.argv[1:])