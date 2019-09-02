from __future__ import division
import random
import futils
import numpy as np

def compose_ep(y_label):
    a = []
    for i, elem in enumerate(y_label):
        a.append(float(input('enter absolute value of element '+str(elem)+' '+str(i+1)+'/'+str(len(y_label))+' in index> ')))
    b = []
    for q in a:
        b.append(q/sum(a))
    return np.matrix(b).reshape(len(b), 1)

def get_ep(y_label, end_points):
    if end_points == 'auto':
        ep = []
        for i in range(len(y_label)):
            ep.append(float(1))
        ep = np.matrix(ep).reshape(len(ep), 1)
    else:
        print(uid)
        ep = compose_ep(y_label)
    return ep


def reformMatrix(row, col, array):
    return np.matrix(array).reshape(row, col)

def svd(S1, S2):
    #########################################
    #Solving n-combinatorial weight function#
    #S1*w = S2                              #
    #w = inv(S1)*S2                         #
    #w^ = max(w)                            #
    #########################################
    try:
        return (S1.I)*S2
    except:
        try:
            #Adding regularization
            row = S1.shape[0]
            col = S1.shape[1]
            J = []
            for i in range(row):
                for k in range(col):
                    J.append(0.05 if i == k else 0.0)
            J = np.matrix(J).reshape(row, col)
            return ((S1+J).I)*S2
        except:
            print("Failed to compute weights matrix...")
            pass

def psvd(IDX, end_points):
    print("Solving matrix...")
    UDX = {}
    #Sequentially solving all uids through all combinations
    for uid in IDX:
        UDX[uid] = {}
        for form in IDX[uid]:
            UDX[uid][form] = {}
            mx = reformMatrix(len(IDX[uid][form][0]), len(IDX[uid][form][1]), IDX[uid][form][2])
            O = get_ep(IDX[uid][form][0], end_points)
            w = svd(mx, O)
            if type(w) != None:
                for i in range(len(IDX[uid][form][1])):
                    UDX[uid][form][IDX[uid][form][1][i]] = float(w[i])
    return UDX
                                                      
