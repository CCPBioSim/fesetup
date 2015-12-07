#
# cheapmap: a simple morph pair mapper
#
# (C) 2015 Hannes H Loeffler, STFC Daresbury
#

from __future__ import print_function

import os

import rdkit.Chem as rd



#MST_NP_FILE = 'scipy_mst.npz'
DOT_FILE = 'mst.dot'
GPICKLE_FILE = 'nx_mst.pickle'
MST_PICKLE_FILE = 'mst.pickle'


# NOTE: the more similar, the smaller the weight must be!
#       0 (or Inf or NaN) means no egde for dense(!) graphs

def tanimoto_score(mol1, mol2):
    """Compute the similarity via Tanimoto fingerprints for mol1 and mol2."""

    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs

    fp1 = FingerprintMols.FingerprintMol(mol1)
    fp2 = FingerprintMols.FingerprintMol(mol2)

    return 1.0 / (DataStructs.FingerprintSimilarity(fp1, fp2) + 1e-15)

def maccs_score(mol1, mol2):
    """Compute the similarity via MACCS fingerprints for mol1 and mol2."""

    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs

    fp1 = rd.MACCSkeys.GenMACCSKeys(mol1)
    fp2 = rd.MACCSkeys.GenMACCSKeys(mol2)

    return 1.0 / (DataStructs.FingerprintSimilarity(fp1, fp2) + 1e-15)


def mcs_score(mol1, mol2):
    """Compute the similarity via the MCS for mol1 and mol2."""

    import math

    from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare

    _fmcs_params = dict(maximizeBonds=False, threshold=1.0, timeout=60,
                        verbose=False, matchValences=False,
                        ringMatchesRingOnly=True, completeRingsOnly=True,
                        atomCompare=AtomCompare.CompareAny,
                        bondCompare=BondCompare.CompareAny)

    mcs = FindMCS( (mol1, mol2), **_fmcs_params)
    pattern = rd.MolFromSmarts(mcs.smartsString)

    # FIXME: deal with multiple matches?
    match1 = mol1.GetSubstructMatch(pattern)

    NA = mol1.GetNumAtoms()
    NB = mol2.GetNumAtoms()
    NMCS = 2 * len(match1)

    # LOMAP, for heavy atoms only
    #beta = 1.0
    #score = math.exp(beta * (NA + NB - NMCS) )

    # simple linear
    score = NA + NB - NMCS + 1

    # Brint&Willet
    #NbA = mol2.GetNumBonds()
    #NbA = mol1.GetNumBonds()
    #mol3 = rd.EditableMol(mol1)

    #for k in range(NA-1, -1, -1):
    #    if k not in match1:
    #        mol3.RemoveAtom(k)

    #NbMCS = mol3.GetMol().GetNumAtoms()

    #score = (float(NMCS + NbMCS) /
    #         float( (NA + NB) * (NbA + NbA) ) )

    return score


valid_methods = {'tanimoto' : tanimoto_score,
                 'maccs' : maccs_score,
                 'mcs' : mcs_score}

def draw_graph(mst, mst_a, mol_names, dir_names, method):

    import networkx as nx

    G = nx.from_scipy_sparse_matrix(mst)

    if method == 'mcs':
        corr = 1
    else:
        corr = 0

    for i, j in zip(mst.nonzero()[0], mst.nonzero()[1]):
        G.edge[i][j]['label'] = '%.1f' % (mst_a[i][j] - corr)
        G.edge[i][j]['len'] = '3.0'

    for n in G.nodes():
        G.node[n]['shape'] = 'box'
        G.node[n]['label'] = ('<'
        '<table border="0" cellspacing="-20" cellborder="0">'
        '<tr><td><img src="%s"/></td></tr>'
        '<tr><td bgcolor="#F0F0F0">%s</td></tr>'
        '</table>>' % (os.path.join(dir_names[n],
                                    mol_names[n] + os.extsep + 'svg'),
                       mol_names[n]) )

    print('Writing networkx graph pickle file %s...' % GPICKLE_FILE)
    nx.write_gpickle(G, GPICKLE_FILE)

    print('Writing DOT file %s...' % DOT_FILE)
    nx.write_dot(G, DOT_FILE)


_mol_params = dict(sanitize=False, removeHs=False)

# FIXME: guard against low scores
#        disallow change in total charge
def calc_MST(filenames, method, do_draw=True, parallel=False):

    from functools import partial

    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree

    import rdkit.Chem.AllChem as ac

    score = valid_methods[method]

    N = len(filenames)
    M = N * (N - 1) / 2
    npout = (M + (100 - M % 100) ) / 100
    simmat = np.zeros(shape=(N,N), dtype=np.float32)

    mols = []
    mol_names = []
    dir_names = []

    print('Reading input files...')

    for filename in filenames:
        mol = rd.MolFromMol2File(filename, **_mol_params)

        dirname = os.path.dirname(filename)
        basename = os.path.splitext(os.path.basename(filename))[0]
        outname = os.path.join(dirname, basename + os.extsep + 'svg')

        mols.append(mol)
        mol_names.append(basename)
        dir_names.append(dirname)

        tmp = ac.Compute2DCoords(mol)

        if do_draw:
            import rdkit.Chem.Draw as draw
            draw.MolToFile(mol, outname, wedgeBonds=False, size=(150,150),
                           fitImage=True, kekulize=False)

    print('Computing similarity matrix using %s...' % method)

    if parallel:
        pool = mp.Pool(mp.cpu_count() )
        map_func = pool.imap
    else:
        map_func = map

    results = []

    for i in range(N-1):
        print('%s...' % mol_names[i])

        partial_func = partial(score, mols[i])
        results.append(map_func(partial_func, mols[i+1:N]) )

    for i, row in enumerate(results):
        simmat[i][i+1:N] = [s for s in row]

    if parallel:
        pool.close()
        pool.join()

    print('similarity score matrix:\n', simmat)

    # NOTE: this removes edges with the larger weight
    mst = minimum_spanning_tree(csr_matrix(simmat))

    cnt = 0
    mst_a = mst.toarray()

    print('\nminimal spanning tree (MST):\n', mst_a)

    print('\nsuggested mappings from MST:')

    for i, j in zip(mst.nonzero()[0], mst.nonzero()[1]):
        cnt += 1
        n1 = mol_names[i]
        n2 = mol_names[j]
        score = mst_a[i][j]

        print('%6i) %s <> %s (%f)\n' % (cnt, n1, n2, score), end='')

    with open(MST_PICKLE_FILE, 'wb') as pfile:
        pickle.dump(mst, pfile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(mol_names, pfile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(dir_names, pfile, pickle.HIGHEST_PROTOCOL)

    return mst, mst_a, mol_names, dir_names



if __name__ == '__main__':

    import argparse
    import sys
    import glob
    import cPickle as pickle

    parser = argparse.ArgumentParser(
        description='Compute the minimal spanning tree (MST) from a set of '
        'mol2 files. The similarity matrix is computed for all the molecules '
        'with a chosen similarity method.  This matrix is then reduced to a '
        'MST path such that all N(N-1)/2 pairs are reduced to a path of N-1 '
        'pairs by minimising the sum of similarity scores.')
    parser.add_argument('mol2_dir', nargs=1,
                        help='directory containing mol2 files '
                        '(only first molecule read)')
    parser.add_argument('-m', '--method', type=str, default=['mcs'], nargs=1,
                        choices=['mcs', 'tanimoto', 'maccs'],
                        help='similarity method used for each molecule pair')
    parser.add_argument('-d', '--draw', action='store_true',
                        help='draw the resulting MST graph (requires PIL, '
                        'pygraphviz)')
    parser.add_argument('-g', '--graph', default=None, nargs=1,
                        help='read graph from previous run and draw it'
                        'pygraphviz)')
    parser.add_argument('-p', '--parallel', action='store_true',
                        help='enable the multiprocessing feature')
    parser.add_argument('--version', action='version', version='%(prog)s 0.2.0')
    parser.add_argument('--tracebacklimit', type=int, default=0, nargs=1,
                        metavar='N',
                        help='set the Python traceback limit to N '
                        '(for debugging)')

    args = parser.parse_args()

    sys.tracebacklimit = args.tracebacklimit


    if args.graph:
        import rdkit.Chem.Draw as draw
        import rdkit.Chem.AllChem as ac

        # FIXME: ugly hack
        for filename in glob.glob('%s/*.mol2' % args.mol2_dir[0]):
            mol = rd.MolFromMol2File(filename, **_mol_params)

            dirname = os.path.dirname(filename)
            basename = os.path.splitext(os.path.basename(filename))[0]
            outname = os.path.join(dirname, basename + os.extsep + 'svg')

            tmp = ac.Compute2DCoords(mol)

            draw.MolToFile(mol, outname, wedgeBonds=False, size=(150,150),
                           fitImage=True, kekulize=False)

        with open(args.graph[0], 'rb') as pfile:
            mst = pickle.load(pfile)
            mol_names = pickle.load(pfile)
            dir_names = pickle.load(pfile)

        mst_a = mst.toarray()
        draw_graph(mst, mst_a, mol_names, dir_names, args.method[0])
    else:
        # FIXME: other file types
        mol2_files = glob.glob('%s/*.mol2' % args.mol2_dir[0])

        if not mol2_files:
            raise IOError('directory %s non-existent or empty' % args.mol2_dir[0])
    
        method = args.method[0]

        if args.parallel:
            import multiprocessing as mp
            print('Running on %i processors...' % mp.cpu_count() )

        mst, mst_a, mol_names, dir_names = calc_MST(mol2_files, method, args.draw,
                                                args.parallel)

        if args.draw:
            draw_graph(mst, mst_a, mol_names, dir_names, method)
