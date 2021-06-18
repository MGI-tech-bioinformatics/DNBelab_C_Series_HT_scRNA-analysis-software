#!/usr/bin/env python3


########
#
#  get.similarityOfBeads.py
#
#  Author: YU Hao (yuhao@genomics.cn)
#          BAI Yinqi (baiyinqi@genomics.cn)
#          LIU weiqing (liuweiqing@genomics.cn)
#
#  Date:   2021-01-20 (v1.0)
#          2021-01-21 (v1.1) To fix up a some of bug
#          2021-03-16 (v1.2) To add a function for filtering the redundancy of pairing-events, 
#                            such as A->B and B->A, A->C and C->A
#          2021-03-17 (v1.3) To modify the extreme limit of m280 count
#          2021-03-18 (v1.4) To add a limit for CPU using
#
########




from datatable import dt, f, g, by, sort
from pandas import Series
from numpy import square, sqrt, sum, dot, percentile
from itertools import product
from sys import argv




########
dt.options.nthreads = 4




########
class WhitelistHash:
    def __init__(self, whitelistFH):
        self._whitelistFH = whitelistFH
        self.whitelistBox = {}


    def hash_whitelist(self):
        for ln in self._whitelistFH:
            ln = ln.strip()

            if 'N' in ln:
                continue

            for i in range(len(ln)):
                new_ln1 = ''.join((ln[:i], 'A', ln[i + 1:]))
                new_ln2 = ''.join((ln[:i], 'C', ln[i + 1:]))
                new_ln3 = ''.join((ln[:i], 'G', ln[i + 1:]))
                new_ln4 = ''.join((ln[:i], 'T', ln[i + 1:]))

                try:
                    self.whitelistBox[new_ln1]
                except KeyError:
                    self.whitelistBox.update({new_ln1: ln})

                try:
                    self.whitelistBox[new_ln2]
                except KeyError:
                    self.whitelistBox.update({new_ln2: ln})

                try:
                    self.whitelistBox[new_ln3]
                except KeyError:
                    self.whitelistBox.update({new_ln3: ln})

                try:
                    self.whitelistBox[new_ln4]
                except KeyError:
                    self.whitelistBox.update({new_ln4: ln})
########




########
class SimilarityCalculation:
    def __init__(self):
        self.similarity = 0.0


    def calCosSimilarity(self, boxX, boxY):
        _Xk = boxX.keys()
        _Yk = boxY.keys()

        _Xv = Series(boxX.values()).rank().to_numpy()
        _Yv = Series(boxY.values()).rank().to_numpy()

        _boxX = dict(zip(_Xk, _Xv))
        _boxY = dict(zip(_Yk, _Yv))

        _Union = _Xk | _Yk
        _box = dict(zip(sorted(_Union), [0] * len(_Union)))

        _boxX_new = _box.copy()
        _boxY_new = _box.copy()

        _boxX_new.update(_boxX)
        _boxY_new.update(_boxY)

        _Xv_new = list(_boxX_new.values())
        _Yv_new = list(_boxY_new.values())

        _Xnorm = _Xv_new / sqrt(sum(square(_Xv_new)))
        _Ynorm = _Yv_new / sqrt(sum(square(_Yv_new)))

        self.similarity = dot(_Xnorm, _Ynorm)
########




sample = argv[1]


#weiqing df = dt.fread('DATA/' + sample + '/CB_UB_count_T1.all.txt', header=False)
#weiqing umi = dt.fread('DATA/' + sample + '/' + sample + '_beads_counts.tsv', header=True)
#weiqing wl_m280 = dt.fread('MISC/m280.barcodes.T.lst', header=False)


df = dt.fread(argv[2], header=False)
umi = dt.fread(argv[3], header=False)
wl_m280 = dt.fread(argv[4], header=False)


droplet_set = set(umi['C0'].to_list()[0])
wl_m280 = wl_m280[:, f.C0].to_list()[0]




wl_m280_1 = set([i[0:10] for i in wl_m280])
wl_m280_2 = set([i[10:20] for i in wl_m280])

wl_m280_1 = WhitelistHash(wl_m280_1)
wl_m280_2 = WhitelistHash(wl_m280_2)

wl_m280_1.hash_whitelist()
wl_m280_2.hash_whitelist()

corrected_m280_1 = [wl_m280_1.whitelistBox[i[0:10]] if i[0:10] in wl_m280_1.whitelistBox else '-NA-' for i in df[:, f.C1].to_list()[0]]
corrected_m280_2 = [wl_m280_2.whitelistBox[i[10:20]] if i[10:20] in wl_m280_2.whitelistBox else '-NA-' for i in df[:, f.C1].to_list()[0]]

corrected_m280 = [''.join(i) if '-NA-' not in i else '-NA-' for i in zip(corrected_m280_1, corrected_m280_2)]


corrected_df = df

corrected_df[:, f.C1] = dt.Frame(corrected_m280)
corrected_df = corrected_df[(f.C1 != '-NA-'), :]
corrected_df = corrected_df[:, dt.sum(f.C0), by([f.C2, f.C1])]




corrected_df_no1 = corrected_df[f.C0 > 1, :]

filtering_lst = corrected_df_no1[:, dt.count(f.C2), by(f.C2)][f.C3 > 1, :]
filtering_lst.key = 'C2'

corrected_df_no1_no1 = corrected_df_no1[:, (f.C2, f.C1, f.C0, g.C3), dt.join(filtering_lst)][f.C3 != None, :]




bead_m280_type_no1_no1 = corrected_df_no1_no1[:, dt.count(f.C2), by(f.C2)]



reclaimed_threshold = 90
extreme_limit = int(percentile(corrected_df_no1_no1['C0'].to_numpy().T[0], reclaimed_threshold))

corrected_df_no1_no1[f.C0 > extreme_limit, f.C0] = extreme_limit




bead_similarity_list = []

droplet_box = {i[0]: dict(zip(list(i[1]['C1']), list(i[1]['C0']))) for i in corrected_df_no1_no1.to_pandas().groupby('C2') if i[0] in droplet_set}

R = 10000

_ = SimilarityCalculation()

for n in range(int(bead_m280_type_no1_no1.nrows / R) + 1):
    bead_list = bead_m280_type_no1_no1[n * R:n * R + R, :]
    bead_list.key='C2'
    corrected_df_no1_no1_subset = corrected_df_no1_no1[:, (f.C2, f.C1, f.C0, g.C3), dt.join(bead_list)][f.C3 != None, :]

    box_bead_m280 = {i[0]: dict(zip(list(i[1]['C1']), list(i[1]['C0']))) for i in corrected_df_no1_no1_subset.to_pandas().groupby('C2')}

    P = product(droplet_box.keys(), box_bead_m280.keys())

    bead_pairList_filtered = [i for i in P if len(droplet_box[i[0]].keys() & box_bead_m280[i[1]].keys()) >= 2]

    for i in bead_pairList_filtered:
#        print(droplet_box[i[0]], box_bead_m280[i[1]])
        _.calCosSimilarity(droplet_box[i[0]], box_bead_m280[i[1]])
        bead_similarity_list.append((i[0], i[1], _.similarity))

bead_similarity_list = dt.Frame(bead_similarity_list)[f.C0 != f.C1, :]




T = 0.8 if len(argv) == 3 and argv[2] == '--atac' else 1.0


box = {}

for i in bead_similarity_list.to_numpy():
    if i[0] not in box:
        box.update({i[0]: {}})

    if i[1] not in box[i[0]]:
        box[i[0]].update({i[1]: i[2]})

    else:
        print('Duplicated value.')


bead_similarity_list_filtered = []
bead_similarity_list_filtered_thread = []


for i in box.keys():
    for j in sorted(box[i].items(), key=lambda x: x[1], reverse=True):
        if j[1] >  T and j[0] not in droplet_set:
            continue
            
        if j[1] <= T and j[0] not in droplet_set:
            break

        bead_similarity_list_filtered.append((i, j[0], j[1]))

bead_similarity_list_filtered = dt.Frame(bead_similarity_list_filtered)





#weiqing bead_similarity_list.to_csv('DATA/' + sample + '/Similarity.all.csv', header=False)
#weiqing bead_similarity_list_filtered.to_csv('DATA/' + sample + '/Similarity.droplet.csv', header=False)


bead_similarity_list.to_csv(argv[5], header=False)
bead_similarity_list_filtered.to_csv(argv[6], header=False)

box = {}

#weiqing with open('DATA/' + sample + '/Similarity.droplet.csv', 'rt') as IN: 
with open(argv[6], 'rt') as IN: 
    for i in IN:
        i = i.strip()
        i = i.split(',')
                
        if tuple([i[0], i[1]]) not in box and tuple([i[1], i[0]]) not in box:
            box.update({tuple([i[0], i[1]]): [i[2]]})

#weiqing with open('DATA/' + sample + '/Similarity.droplet.filtered.csv', 'wt') as OU:
with open(argv[7], 'wt') as OU:
    for i in box:
        n = list(i)
        n.extend(box[i])
        print(','.join(n), file=OU)
