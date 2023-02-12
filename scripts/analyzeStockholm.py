#!/usr/bin/env python3
'''Analyze correlation of basepairs from an annotated stockholm alignment file'''
import re, sys, csv, pprint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import AlignIO
from collections import Counter
from itertools import permutations
from Bio.Align import MultipleSeqAlignment

def getBasePairs(seq, bpPos):
    for bp in bpPos:
        yield (seq[bp[0]], seq[bp[1]])

def is_base_pair(x, y):
    x, y = x.lower(), y.lower()
    return (x=='a' and y=='u' or
            x=='u' and y=='a' or
            x=='c' and y=='g' or
            x=='g' and y=='c' or
            x=='g' and y=='u' or
            x=='u' and y=='g' )

def organizeBPsBYindex(align: MultipleSeqAlignment, indexLoc:tuple):
    '''Given an alignment  and start-end location as a tuple for organizing,
    returns the base-pairs from bpPos organized by the sequence of indexLoc'''
    aLoop2554_2555_to_BPs = dict()

    openPos = [i for i, x in enumerate(align.column_annotations["secondary_structure"]) if x == '(']
    closePos = [i for i, x in enumerate(align.column_annotations["secondary_structure"]) if x == ')']
    bpPos = list(zip(openPos, reversed(closePos)))
    for seqEntry in align:
        if re.match(r'^LSU.*', seqEntry.id):
            continue
        nucl_2554_2555_str = str(seqEntry.seq[indexLoc[0]:indexLoc[1]])
        if nucl_2554_2555_str not in aLoop2554_2555_to_BPs.keys():
            aLoop2554_2555_to_BPs[nucl_2554_2555_str] = []
        basePairLetters = list(getBasePairs(seqEntry.seq, bpPos))
        aLoop2554_2555_to_BPs[nucl_2554_2555_str].append((seqEntry.id, seqEntry.description, basePairLetters))
    return aLoop2554_2555_to_BPs

def organizeBPsByPhyloAndALoopID(bpPosition, aLoop2554_2555_to_BPs, possibleRNAcombinations, phyloType='phylum'):
    '''Organize BPs by phylogeny and by A-loop identity
    bpPosition: int, position of base-pair in alignment
    aLoop2554_2555_to_BPs: dict, keys are A-loop identities, values are lists of tuples (phylo, seqID, basePairLetters)
    possibleRNAcombinations: list, list of possible RNA combinations
    phyloType: str, type of phylogeny to use, either "phylum" or "class"'''
    phyloTranslator = {'phylym': 2, 'class': 3}
    if phyloType not in phyloTranslator.keys():
        raise ValueError(f'phyloType must be either "phylum" or "class"')
    
    numberOfBasePairTypes, bpTypeToALoopIdAndPhylogeny, aLoopToPhylo, bpTypeToPhylo, phylos  = dict(), dict(), dict(), dict(), set()
    for k, v in aLoop2554_2555_to_BPs.items():
        for item in v:
            splitPhylo = item[1].split(';')
            if len(splitPhylo) > phyloTranslator[phyloType]-1:
                phylos.add(splitPhylo[phyloTranslator[phyloType]-1])
            else:
                phylos.add(splitPhylo[-1])
        bpTypes = Counter([f'{x[2][bpPosition][0]}{x[2][bpPosition][1]}' for x in v])
        specTypes = dict(Counter([f'{";".join(x[1].split(";")[0:phyloTranslator[phyloType]])}%${x[2][bpPosition][0]}{x[2][bpPosition][1]}' for x in v]))
        aLoopToPhylo[k] = dict(Counter([f'{";".join(x[1].split(";")[0:phyloTranslator[phyloType]])}' for x in v]))
        bpTypesOut = list()
        for bp in possibleRNAcombinations:
            if bp in bpTypes:
                bpTypesOut.append(bpTypes[bp])
            else:
                bpTypesOut.append(0)
        numberOfBasePairTypes[k] = bpTypesOut
        for phyloBP, number in specTypes.items():
            phylo, bp = phyloBP.split('%$')
            if bp not in bpTypeToALoopIdAndPhylogeny.keys():
                bpTypeToALoopIdAndPhylogeny[bp] = dict()
            if k not in bpTypeToALoopIdAndPhylogeny[bp].keys():
                bpTypeToALoopIdAndPhylogeny[bp][k] = list()
            bpTypeToALoopIdAndPhylogeny[bp][k].append((phylo,number))
            if bp not in bpTypeToPhylo.keys():            
                bpTypeToPhylo[bp] = dict()
            if phylo not in bpTypeToPhylo[bp].keys():
                bpTypeToPhylo[bp][phylo] = 0
            bpTypeToPhylo[bp][phylo] += number
    return numberOfBasePairTypes,bpTypeToALoopIdAndPhylogeny,aLoopToPhylo,bpTypeToPhylo,phylos

def plotPieFromListOfTups(listWithTups, axisObj, colorDict=None, splitBy=None):
    colors = None
    if splitBy:  
        namesPh = [x[0].split(splitBy)[-1] for x in listWithTups]
    else:
        namesPh = [x[0] for x in listWithTups]
    if colorDict:
        colors = [colorDict[name] for name in namesPh]
    vals = [x[1] for x in listWithTups]
    axisObj.pie(vals, colors=colors)
    return axisObj

def drawPhyloDistributionPies(bpPosition, bpTypeToALoopIdAndPhylogeny, aLoopToPhylo, bpTypeToPhylo, df, phylos, tempCSV):

    matchesCapital = re.findall(r'(?:(?<=^)|(?<=[^.]))\s+([A-Z][a-z]+)', ' '.join(sorted(phylos)))
    matchesNonCapital = re.findall(r'(?:(?<=^)|(?<=[^.]))\s+([a-z]+)', ' '.join(sorted(phylos)))
    matchesCapital.remove('Candidatus')
    matchesNonCapital.append('Candidatus')
    matchesNonCapital.append('1286008:1289191')
    color_range = list(np.linspace(0, 1, len(matchesCapital), endpoint=False))

    #Integrate colors here based on temperature
    for name in matchesCapital:
        if name in tempCSV.keys():
            pass

    colors = [plt.cm.tab20.colors,plt.cm.tab20b.colors]
    flattenedColors = [item for sublist in colors for item in sublist]

    colors = [flattenedColors[i] for i,x in enumerate(matchesCapital)]
    color_dict = dict(zip(matchesCapital, colors))    
    for x in matchesNonCapital:
        color_dict[x] = (0, 0, 0, 1)
    fig, axs = plt.subplots(nrows=len(df.index)+1, ncols=len(df.columns)+1)
    for rowI, RowAx in enumerate(axs):
        for colI, ax in enumerate(RowAx):
            if colI >= len(df.columns) and rowI >= len(df.index):
                ax.axis('off')
                continue
            if colI >= len(df.columns):
                #do column relevant ones
                plotPieFromListOfTups([(k, v) for k, v in aLoopToPhylo[df.index[rowI]].items()], ax, splitBy=';', colorDict=color_dict)
                continue
            if rowI >= len(df.index):
                #do the rows
                plotPieFromListOfTups([(k, v) for k, v in bpTypeToPhylo[df.columns[colI]].items()], ax, splitBy=';', colorDict=color_dict)
                continue
            if df.columns[colI] not in bpTypeToALoopIdAndPhylogeny.keys():
                ax.axis('off')
                continue
            if df.index[rowI] not in bpTypeToALoopIdAndPhylogeny[df.columns[colI]].keys():
                ax.axis('off')
                continue
            axDone = plotPieFromListOfTups(bpTypeToALoopIdAndPhylogeny[df.columns[colI]][df.index[rowI]], ax, splitBy=';', colorDict=color_dict)
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in color_dict.values()]
    plt.legend(markers, color_dict.keys(), numpoints=1, loc=(1.04,0))
    plt.savefig(f'piesBPpos{bpPosition+1}.svg', bbox_inches="tight")

def plotHeatMapFromDF(df, bpPosition, orientation='vertical'):
    heatmap = plt.imshow(df)
    plt.xticks(range(len(df.columns.values)), df.columns.values)
    plt.yticks(range(len(df.index)), df.index)
    plt.xlabel(f'Base-pair at position {bpPosition+1}')
    plt.ylabel('2554 2555 identity')
    cbar = plt.colorbar(mappable=heatmap, orientation=orientation)  

    for i, ikey in enumerate(df.columns):
        for j, jkey in enumerate(df.index):
            text = plt.text(i, j, df[ikey][jkey],
                           ha="center", va="center", color="w")
    

    plt.savefig(f'bpPos{bpPosition+1}.svg')
    return True

def main():
    align = AlignIO.read("alns/H92_fullArch.sto", "stockholm")
    alignWithDescriptions = AlignIO.read("alns/H92_fullArch.fa", "fasta")
    
    temperatureCSV = dict()
    with open("200617_TEMPURA.csv", 'r') as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            #Name, taxID, SuperK, Phylum, Class, Tmin, Topt_Ave, Tmax, Tmax_Tmin
            if line[6] == '':
                line[6] = line[0]
            if line[5] not in temperatureCSV.keys():
                temperatureCSV[line[5]] = []
            temperatureCSV[line[5]].append(float(line[15]))
            #temperatureCSV.append([line[0].replace('"',''), line[1], line[3], line[4], line[5], line[14], line[15], line[18], line[19]])

    bpPosition = int(sys.argv[1])-1
    for i, seq in enumerate(alignWithDescriptions):
        if seq.id == align[i].id and align[i].description != seq.description:
            align[i].description = seq.description.split()[1]
            
    aLoop2554_2555_to_BPs = organizeBPsBYindex(align, (7,9))
    possibleRNAcombinations = [''.join(x) for x in list(permutations(['A','U','C','G'], 2))]    
    
    numberOfBasePairTypes, bpTypeToALoopIdAndPhylogeny, aLoopToPhylo, bpTypeToPhylo, phylos = organizeBPsByPhyloAndALoopID(bpPosition, aLoop2554_2555_to_BPs, possibleRNAcombinations, 'class')


    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(bpTypeToALoopIdAndPhylogeny)
    df = pd.DataFrame(numberOfBasePairTypes)
    df = df[df.sum(axis=1) > 5]
    df = df.transpose()
    df.columns = list(np.array(possibleRNAcombinations)[df.columns])
    df = df[df.sum(axis=1) > 5]

    # with open('../ggkbase_color_scheme_extended.csv', mode='r') as inp:
    #     reader = csv.reader(inp)
    #     phColors = {rows[0]:rows[1] for rows in reader}
    drawPhyloDistributionPies(bpPosition, bpTypeToALoopIdAndPhylogeny, aLoopToPhylo, bpTypeToPhylo, df, phylos, temperatureCSV)
    
    plt.cla()
    plt.clf()

    plotHeatMapFromDF(df, bpPosition)

if __name__ == '__main__':
    sys.exit(main())
