#!/usr/bin/env python3
'''Given a SILVA .tre and map files from https://www.arb-silva.de/no_cache/download/archive/release_138.1/Exports/taxonomy/
as well as TEMPURA csv (from http://togodb.org/db/tempura), generates a color file for ITOL and adjusted tree with names and if needed truncations.'''
import csv, sys, argparse
from ete3 import NCBITaxa, Tree
from Bio import AlignIO
import matplotlib.cm as cm
from matplotlib.colors import Normalize, rgb2hex

from analyzeStockholm import organizeBPsBYindex

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

def main(commandline_args):
    comm_args = create_and_parse_argument_options(commandline_args)
    tree = Tree(comm_args.treeFile, format=1)
    align = AlignIO.read(comm_args.alignmentFile, "stockholm")
    nuclOne, nuclTwo = int(comm_args.nucleotideIndex[0]), int(comm_args.nucleotideIndex[1])

    nameMapper, taxonomyMapper = dict(), dict()
    with open(comm_args.mapFile) as mapFile:                                                                                          
        reader = csv.DictReader(mapFile, delimiter='\t')
        next(reader)
        for entry in reader:
            if entry['taxid'] not in nameMapper.keys():
                nameMapper[entry['taxid']] = list()
            nameMapper[entry['taxid']].append((entry['organism_name'], entry['primaryAccession'], entry['path'])) 
    
    with open(comm_args.treeFile.replace('.tre','.map')) as mapFile:                                                                                          
        reader = csv.reader(mapFile, delimiter='\t')
        for entry in reader:
            taxonomyMapper[entry[0]] = entry[1]
    
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
                    
                    #Be careful with this line, it is a bit of a hack to get the nucleotide sequence from the alignment
                    #it won't work for proper ranges of sequences, but it will work for the current case
                    if len(node.children) == 0:
                        for sequenceEntry in nameMapper[node.name]:                      
                            node.add_child(name=sequenceEntry[0],  dist=0)
                    accessionToNucl[accession].append(str(align[index].seq[nuclOne:nuclTwo]))
                    node.nuclOfInterest.append(str(align[index].seq[nuclOne:nuclTwo]))
            node.add_features(silvaID=node.name)
            #The name from .map
            #Special one for all the matches in nameMapper
            node.name = taxonomyMapper[node.name]
        else:
            node.delete()

    if comm_args.taxonomyName:
        nodeOfInterest = tree.search_nodes(name=comm_args.taxonomyName)[0]
        truncatedTree = nodeOfInterest.detach()
    else:
        truncatedTree = tree   

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
    norm = Normalize(vmin=-50, vmax=max(toptList), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.turbo)

    with open(f'./temperatureTrees/temperatureLabelColors_{comm_args.taxonomyName}.txt', 'w') as f:
        f.write('DATASET_STYLE\r\nSEPARATOR COMMA\r\nDATASET_LABEL,SILVA tree\r\nCOLOR,#ffff00\r\nDATA\r\n')
        for node in truncatedTree.traverse():
            labelColor = '#ffffff'
            if 'topt_ave' in node.features:
                labelColor = rgb2hex(mapper.to_rgba(node.topt_ave))
            f.write(f'{node.name},label,node,#000000,1,normal,{labelColor}\r\n')

    truncatedTree.write(format=1, outfile=f'./temperatureTrees/truncatedTree_{comm_args.taxonomyName}.nwk')

    #Get the BP data
    
    alignWithDescriptions = AlignIO.read(comm_args.alignmentFile.replace(".sto",".fa"), "fasta")


    for i, seq in enumerate(alignWithDescriptions):
        if seq.id == align[i].id and align[i].description != seq.description:
            align[i].description = seq.description.split()[1]
    bpData = organizeBPsBYindex(align, (int(comm_args.nucleotideIndex[0]), int(comm_args.nucleotideIndex[1])))
    pass


if __name__ == "__main__":
    main(sys.argv[1:])