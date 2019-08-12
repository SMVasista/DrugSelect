from __future__ import division
import futils
import numpy as np

def compose_ep(y_label):
    a = []
    for elem in y_label:
        a.append(float(input('enter absolute value of element'+str(elem)+' in index> ')))
    b = []
    for q in a:
        b.append(q/sum(a))
    return np.matrix(b).reshape(len(b), 1)

def psvd(IDX, end_points):
    UDX = {}
    #Sequentially solving all uids through all combinations
    for uid in IDX:
        UDX[uid] = {}
        for form in IDX[uid]:
            #Solving n-combinatorial weight function
            #U*w = I
            #w = inv(U)*I
            #w^ = max(w)
            mx = np.matrix(IDX[uid][form][2]).reshape(len(IDX[uid][form][0]), len(IDX[uid][form][1]))
            if end_points == 'auto':
                ep = []
                for i in range(len(IDX[uid][form][0])):
                    ep.append(float(1))
                ep = np.matrix(ep).reshape(len(ep), 1)
            else:
                ep = compose_ep(IDX[uid][form][0])
            #######################################
            w = mx.I*ep
            #######################################
            for i in range(len(IDX[uid][form][1])):
                UDX[uid][form][IDX[uid][form][1][i]] = w[i]
    return UDX
                                                      
