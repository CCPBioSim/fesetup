#
# cheapmap: a simple morph pair mapper
#
# (C) 2015 Hannes H Loeffler, STFC Daresbury
#

from __future__ import print_function

import os
import sys
import glob
import math

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import rdkit.Chem as rd
import rdkit.Chem.AllChem as ac
import rdkit.Chem.Draw as draw
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare


beta = 1.0

_fmcs_params = dict(maximizeBonds = False, threshold = 1.0, timeout=1,
                    verbose = False, matchValences = False,
                    ringMatchesRingOnly = True, completeRingsOnly = True,
                    atomCompare = AtomCompare.CompareAny,
                    bondCompare = BondCompare.CompareAny)

_mol_params = dict(sanitize = False, removeHs = False)

# FIXME: guard afainst low scores
def calc_MST(filenames, sim_method):

    N = len(filenames)
    M = N * (N - 1) / 2
    npout = (M + (100 - M % 100) ) / 100
    cgraph = np.zeros(shape=(N,N), dtype=np.float32)

    mols = []
    mol_names = []

    print('Reading input files...')

    for filename in filenames:
        mol = rd.MolFromMol2File(filename, **_mol_params)
        mols.append(mol)
        basename = os.path.splitext(os.path.basename(filename))[0]
        mol_names.append(basename)
        outname = basename + os.extsep + 'png'
        tmp = ac.Compute2DCoords(mol)
        draw.MolToFile(mol, outname, wedgeBonds=False,
                       fitImage=True, kekulize=False)

    print('Computing similarity matrix using %s...' % sim_method)

    cnt = 0 

    for i in range(N-1):
        mol1 = mols[i]

        for j in range(i+1, N):
            mol2 = mols[j]

            # FIXME: disallow change in total charge

            if sim_method == 'tanimoto':
                fp1 = FingerprintMols.FingerprintMol(mol1)
                fp2 = FingerprintMols.FingerprintMol(mol2)
                score = DataStructs.FingerprintSimilarity(fp1, fp2)
            elif sim_method == 'maccs':
                fp1 = rd.MACCSkeys.GenMACCSKeys(mol1)
                fp2 = rd.MACCSkeys.GenMACCSKeys(mol2)
                score = DataStructs.FingerprintSimilarity(fp1, fp2)
            elif sim_method == 'mcs':
                mcs = FindMCS( (mol1, mol2), **_fmcs_params)
                pattern = rd.MolFromSmarts(mcs.smartsString)

                # FIXME: deal with multiple matches?
                match1 = mol1.GetSubstructMatch(pattern)

                NA = mol1.GetNumAtoms()
                NB = mol2.GetNumAtoms()
                NMCS = 2 * len(match1)

                # NOTE: the more similar, the smaller the weight must be!
                #       0 (or Inf or NaN) mean no egde for dense(!) graphs
                #
                #       LOMAP, for heavy atoms only
                #score = math.exp(beta * (NA + NB - NMCS) )

                # simple linear
                score = NA + NB - NMCS + 1

                # Brint&Willethttp://effbot.org/imagingbook/imagetk.htm
                #NbA = mol2.GetNumBonds()
                #NbA = mol1.GetNumBonds()
                #mol3 = rd.EditableMol(mol1)
 
                #for k in range(NA-1, -1, -1):
                #    if k not in match1:
                #        mol3.RemoveAtom(k)

                #NbMCS = mol3.GetMol().GetNumAtoms()

                #score = (float(NMCS + NbMCS) /
                #         float( (NA + NB) * (NbA + NbA) ) )


                # NOTE: triangle or symmetric matrix? Possible multiplicity:
                #       either accept this or change scoring function
            else:
                raise ValueError('Unknown similarity method: %s'
                                 % sim_method)

            cgraph[i][j] = score
            cnt += 1

            if not cnt % npout:
                print('%i...' % cnt)

    print('similarity score matrix:\n', cgraph)

    # NOTE: this removes edges with the larger weight
    mst = minimum_spanning_tree(csr_matrix(cgraph) )

    mst_a = mst.toarray()
    print('\nminimal spanning tree (MST):\n', mst_a)

    cnt = 0
    edge_labels = {}
    print('\nsuggested mappings from MST:')

    for i, j in zip(mst.nonzero()[0], mst.nonzero()[1]):
        cnt += 1
        n1 = mol_names[i]
        n2 = mol_names[j]
        score = mst_a[i][j]
        print('%6i) %s <> %s (%f)\n' % (cnt, n1, n2, score), end='')
        edge_labels[i,j] = '%.4f' % score

    # draw picture for each molecule
    #imgs = []

    #for mol in mols:
    #   tmp = ac.Compute2DCoords(mol)
    #   imgs.append(draw.MolToMPL(mol,kekulize=False) )  # matplotlib.figure.Figure

    #img = draw.MolsToGridImage(mols, molsPerRow=20, subImgSize=(150,150),
    #                           legends=mol_names, kekulize=False)
    #img.save('test.png')

    G = nx.from_scipy_sparse_matrix(mst)

    for n in G:
        img = mpimg.imread(mol_names[n] + os.extsep + 'png')
        G.node[n]['image'] = img

    pos = nx.spring_layout(G, scale = 1.0, iterations=2000)
    #print(pos)
    node_labels = {}

    for i, name in enumerate(mol_names):
        node_labels[i] = name

    plt.figure(0, figsize=(N,N), dpi=80)

    #nx.draw_networkx_nodes(G, pos, node_size=20, node_color='r', node_shape='o')
    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_edge_labels(G, pos, edge_labels)
    #nx.draw_networkx_labels(G, pos, node_labels)
    nx.draw_networkx_labels(G, pos)

    ax = plt.gca()
    fig = plt.gcf()
    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform
    imsize = 0.05

    for n in G.nodes():
        (x, y) = pos[n]
        xx, yy = trans((x, y)) # figure coordinates
        xa, ya = trans2((xx, yy)) # axes coordinates

        img =  G.node[n]['image']

        a = plt.axes([xa - imsize / 2.0, ya - imsize / 2.0, imsize, imsize])
        a.imshow(img)
        a.set_aspect('equal')
        a.axis('off')

    plt.axis('off')
    plt.savefig('mst.svg')



if __name__ == '__main__':

    if len(sys.argv) < 3:
        print('Usage: %s similarity[tanimoto|maccs|mcs] dir_with_mol2' %
              sys.argv[0], file=sys.stderr)
        exit(1)

    # FIXME: other file types
    mol2_files = glob.glob('%s/*.mol2' % sys.argv[2])

    calc_MST(mol2_files, sys.argv[1])
