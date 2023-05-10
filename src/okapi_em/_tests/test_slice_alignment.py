from okapi_em import slice_alignment
import numpy as np
from numpy.random import randint

def test_cl_affineTransform():

    #Create linear transform matrix and a translation matrix

    points = randint(0,10,size=2) #X,Y
    print(points)
    # identity style transfrom
    lin_tr = np.eye(2)
    transl_tr=np.array([0,0])

    myAf0 = slice_alignment.affineTransform(lin_tr, transl_tr)

    #Apply transform and check result
    p_res = myAf0.applyToPoint(points)

    assert p_res[0]==points[0]
    assert p_res[1]==points[1]

    #rotation 90 degrees
    lin_tr = np.array([[ 0, -1], [1, 0]])
    myAf1 = slice_alignment.affineTransform(lin_tr, transl_tr)
    p_res = myAf1.applyToPoint(points)
    assert p_res[0]==-points[1]
    assert p_res[1]==points[0]

    #Translation
    transl_tr=randint(10,20,size=2)
    myAf2 = slice_alignment.affineTransform(lin_tr, transl_tr)
    p_res = myAf2.applyToPoint(points)
    assert p_res[0]==-points[1]+transl_tr[0]
    assert p_res[1]==points[0]+transl_tr[1]