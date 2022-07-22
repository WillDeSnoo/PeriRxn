#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 12:01:32 2022

@author: william
"""

import numpy as np

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def rot_m(xyz,rm):
    lrxyz=[]
    for a in xyz:
        lrxyz.append(np.dot(rm,a))
    return np.asarray(lrxyz)

v1=np.array([1,2,3])
v2=np.array([2,1,3])

rm=rotation_matrix_from_vectors(v1, v2)
print(rm)

v3=np.dot(v2,rm)

print(v1)
print(v2)
print(v3)


xyz=np.asarray([[1,2,3],
     [1,0,1],
     [2,1,0]])

b=rot_m(xyz,rm)
print(b)




        