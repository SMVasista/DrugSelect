from __future__ import division
import sys, os, re
import numpy as np

#Sample input data
Pcomp = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8']
C = [['P1','high','escape'], ['P2','high','accelerate'], ['P3','high','escape'], ['P7', 'low', 'escape']]
H1 = [['P1', 'direct', 'inhibit'], ['P4', 'direct', 'activate']]
H2 = [['P4', 'direct', 'activate']]
H3 = [['P2','direct','activate'], ['P3','indirect','inhibit']]
I = {'P1': 0.24, 'P2': 0.13, 'P3': 0.01, 'P4': 0.001, 'P5': 0.55, 'P6': 0.12, 'P7': 0.25, 'P8': 0.11}

#Translated above information as per the previous scoring techniques etc
I_ = np.array([[0.24], [0.13], [0.01], [0.001], [0.55], [0.12], [0.25], [0.11]])
H1 = np.array([[-1], [0], [0], [1], [0], [0], [0], [0]])
H2 = np.array([[0], [0], [0], [1], [0], [0], [0], [0]])
H3 = np.array([[0], [1], [-0.5], [0], [0], [0], [0], [0]])

#Expectation matrix is the conditional complement of C component
C = np.array([[-0.24], [0.07], [-0.005], [0], [0], [0], [0.25], [0]])
#Therefore the expectation matrix will be
E = np.array([[-1], [1], [-1], [0], [0], [0], [1], [0]])

# J = I_ * H(i) #I_ not considered for simplification
# R = 1/(n+2) *(J*E)
#Performing element-wise operation on the components to maximize function R
#1D ranks

print("\n~~~~1D~~~~")
print("H1:", (1/2)*(sum((H1)*E)))
print("H2:", (1/2)*(sum((H2)*E)))
print("H3:", (1/2)*(sum((H3)*E)))

#2D ranks
print("\n~~~~2D~~~~")
print("H1_H2:", (1/3)*sum((H1+H2)*E))
print("H2_H3:", (1/3)*sum((H2+H3)*E))
print("H3_H1:", (1/3)*sum((H3+H1)*E))

#3D ranks
print("\n~~~~3D~~~~")
print("H1_H2_H3:", (1/4)*sum((H1+H2+H3)*E))

K = {'H1': (1/2)*(sum((H1)*E))[0], 'H2': (1/2)*(sum((H2)*E))[0], 'H3': (1/2)*(sum((H3)*E))[0], 'H1_H2': (1/3)*sum((H1+H2)*E)[0], 'H2_H3': (1/3)*sum((H2+H3)*E)[0], 'H3_H1': (1/3)*sum((H3+H1)*E)[0], 'H1_H2_H3': (1/4)*sum((H1+H2+H3)*E)[0]}

#Penalizing off-target effects of Hx
#This can be done in the following manner
#P(Hx) = 5/(5 + (nonzero(Hx) - nonzero(Hx*E))) #syntax 5/(5 + (np.count_nonzero(H) - np.count_nonzero(H*E)))
#and this P(Hx) multiplied into the final rank value

print("\n######################")
print("Results - sorted order")
print("######################\n")
for elem in sorted(K.items(), key=lambda x: x[1], reverse=True):
	print(elem[0], ':', elem[1])

print("\n######################\nBased on our expectation H3 and H3+H1 are the best compliments to the drug C - This is visible from the results also. H3+H1 is still ranked higher because  the  degree of complement is very high in the combination")

