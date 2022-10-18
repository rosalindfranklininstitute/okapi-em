'''
Copyright 2022 Rosalind Franklin Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

# Code originally developed by Chloe Chen
# Adapted here to run on okapi-em napari plugin

import numpy as np
import cv2 as cv #opencv package
from scipy.optimize import least_squares


ALIGNMENT_METHOD_DEFAULT = {
    'translation':True, # (can be overriden by affine setting below)
    'rotation':False, 'shearing_x':True, 'shearing_y':False, # If rotation is true the shearings are ignored
    'scaling':False, 'stretching_x':False, 'stretching_y':True, # If scaling is true, the stretchings are ignored
    'affine':False # If affine is true all the others are ignored. This setting mean that all parameters in the affine matrix are free to adjust in the optimization
}


def featureMatchingDS(img1, img2):
    # Initiate SIFT detector
    sift = cv.SIFT_create()
    # find the keypoints and descriptors with SIFT
    kp1, des1 = sift.detectAndCompute(img1,None)
    kp2, des2 = sift.detectAndCompute(img2,None)
    
    # BFMatcher with default params
    bf = cv.BFMatcher()
    matches = bf.knnMatch(des1,des2,k=2)

    good = []
    good_no_list = []
    for m, n in matches:
        #if m.distance < 0.5*n.distance:
        #if m.distance < 0.75*n.distance:
        if m.distance < 0.65*n.distance:
            good.append([m])
            good_no_list.append(m)

    reference_pts = np.float32([kp1[m.queryIdx].pt for m in good_no_list]).reshape(-1, 1, 2)     
    new_pts = np.float32([kp2[m.trainIdx].pt for m in good_no_list]).reshape(-1, 1, 2)

    x_ref, y_ref = reference_pts[:, 0, 0], reference_pts[:, 0, 1]
    x_new, y_new = new_pts[:, 0, 0], new_pts[:, 0, 1]

    slope = np.arctan2((y_new - y_ref), (x_new - x_ref))
    dist = np.sqrt((y_new - y_ref) ** 2 + (x_new - x_ref) ** 2)
    
    points_ref = [img1[x, y] for x, y in zip(y_ref.astype(int), x_ref.astype(int))]
    points_new = [img2[x, y] for x, y in zip(y_new.astype(int), x_new.astype(int))]
    
    #slope_select = np.logical_and(slope > np.median(slope) - 0.08, slope < np.median(slope) + 0.08)
    dist_select = np.logical_and(dist > np.median(dist) * 0.95, dist < np.median(dist) * 1.05)
    points_select = np.logical_and(np.array(points_ref) > 1, np.array(points_new) > 1)
    
    select = dist_select * points_select

    matched_filtered = np.array(good)[select]
    reference_pts = reference_pts[select]    
    new_pts = new_pts[select]
    
    print(f'{len(matched_filtered)} good matches found')

    # cv.drawMatchesKnn expects list of lists as matches.
    img_match_gauss = cv.drawMatchesKnn(img1,kp1,img2,kp2,matched_filtered,None,flags=cv.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS)
    #plt.figure(figsize = (20,20))
    #plt.imshow(img_match_gauss)
    
    return reference_pts, new_pts

def transform(method, param, points):
    x0 = points[0]
    y0 = points[1]
    
    #rotation or shearing
    if method['rotation'] == True:
        x1 = np.cos(param[2]) * x0 - np.sin(param[2]) * y0
        y1 = np.sin(param[2]) * x0 + np.cos(param[2]) * y0
    
    elif method['shearing_x'] == True:
        x1 = x0 + param[2] * y0
        y1 = y0
        
    elif method['shearing_y'] == True:
        x1 = x0
        y1 = param[2] * x0 + y0
    
    else:
        x1 = x0
        y1 = y0
    
    #scaling or stretching
    if method['scaling'] == True:
        x1 = x1 * param[3]
        y1 = y1 * param[3]
        
    elif method['stretching_x'] == True and method['stretching_y'] == False:
        x1 = x1 * param[3]
        y1 = y1
    
    elif method['stretching_y'] == True and method['stretching_x'] == False:
        x1 = x1
        y1 = y1 * param[3]
        
    elif method['stretching_x'] == True and method['stretching_y'] == True:
        x1 = x1 * param[3]
        y1 = y1 * param[4]
        
        
    #affine 
    if method['affine'] == True:
        x1 = param[2] * x0 + param[3] * y0
        y1 = param[4] * x0 + param[5] * y0
        
        
    #translation    
    if method['translation'] == True:
        x1 += param[0]
        y1 += param[1]
        
    return x1, y1

def findMatrix(method, params):
    if method['affine'] == True:
        T = np.array([[params[2], params[3], params[0]], [params[4], params[5], params[1]]])
        
    elif method['scaling'] == True:
        if method['rotation'] == True:
            T = np.array([[params[3] * np.cos(params[2]), -params[3] * np.sin(params[2]), params[0]], 
                          [params[3] * np.sin(params[2]), params[3] * np.cos(params[2]), params[1]]])
        elif method['shearing_x'] == True:
            T = np.array([[params[3], params[3] * params[2], params[0]], [0, params[3], params[1]]])
        elif method['shearing_y'] == True:
            T = np.array([[params[3], 0, params[0]], [params[3] * params[2], params[3], params[1]]])
        else:
            T = np.array([[params[2], 0, params[0]], [0, params[2], params[1]]])
            
    elif method['stretching_x'] == True:
        if method['rotation'] == True:
            T = np.array([[params[3] * np.cos(params[2]), -np.sin(params[2]), params[0]], 
                          [params[3] * np.sin(params[2]), np.cos(params[2]), params[1]]])
        elif method['shearing_x'] == True:
            T = np.array([[params[3], params[2], params[0]], [0, 1, params[1]]])
        elif method['shearing_y'] == True:
            T = np.array([[params[3], 0, params[0]], [params[3] * params[2], 1, params[1]]])
        else:
            T = np.array([[params[2], 0, params[0]], [0, 1, params[1]]])
            
    elif method['stretching_y'] == True:
        if method['rotation'] == True:
            T = np.array([[np.cos(params[2]), -params[3] * np.sin(params[2]), params[0]], 
                          [np.sin(params[2]), params[3] * np.cos(params[2]), params[1]]])
        elif method['shearing_x'] == True:
            T = np.array([[1, params[3] * params[2], params[0]], [0, params[3], params[1]]])
        elif method['shearing_y'] == True:
            T = np.array([[1, 0, params[0]], [params[2], params[3], params[1]]])
        else:
            T = np.array([[1, 0, params[0]], [0, params[2], params[1]]])
            
    else:
        if method['rotation'] == True:
            T = np.array([[np.cos(params[2]), -np.sin(params[2]), params[0]], [np.sin(params[2]), np.cos(params[2]), params[1]]])
        elif method['shearing_x'] == True:
            T = np.array([[1, params[2], params[0]], [0, 1, params[1]]])
        elif method['shearing_y'] == True:
            T = np.array([[1, 0, params[0]], [params[2], 1, params[1]]])
        else:
            T = np.array([[1, 0, params[0]], [0, 1, params[1]]])
            
    affine = T[:, :2]
    translate = T[:, 2]

    return T, affine, translate

def residual(ref, new, method):
    def r(t):
        rx = transform(method, t, new)[0] - ref[0]
        ry = transform(method, t, new)[1] - ref[1]
        return np.append(rx, ry)
    return r

def findInitialandBounds(method, margin_t=200, margin_a=0.01):
    if method['affine'] == True:
        t0 = np.array([0, 0, 1, 0, 0, 1])
        bounds = ([-margin_t, -margin_t, 1 - margin_a, -margin_a, -margin_a, 1 - margin_a]
                  , [margin_t, margin_t, 1 + margin_a, margin_a, margin_a, 1 + margin_a])
    
    elif method['scaling'] == True or method['stretching_x'] == True or method['stretching_y'] == True:
        if method['rotation'] == True or method['shearing_x'] == True or method['shearing_y'] == True:
            t0 = np.array([0, 0, 0, 1])
            bounds = ([-margin_t, -margin_t, -margin_a, 1 - margin_a], [margin_t, margin_t, margin_a, 1 + margin_a])
        else:
            t0 = np.array([0, 0, 1])
            bounds = ([-margin_t, -margin_t, 1 - margin_a], [margin_t, margin_t, 1 + margin_a])
    
    else:
        if method['rotation'] == True or method['shearing_x'] == True or method['shearing_y'] == True:
            t0 = np.zeros(3)
            bounds = ([-margin_t, -margin_t, -margin_a], [margin_t, margin_t, margin_a])
        else:
            t0 = np.zeros(2)
            bounds = ([-margin_t, -margin_t], [margin_t, margin_t])
            
    return t0, bounds

class affineTransform():
    def __init__(self,lintransf, transltransf):
        self.lintransf = lintransf.copy()
        self.transltransf= transltransf.copy()
    
    @staticmethod
    def from2x3mat(mat):
        lintransf = mat[:, :2].astype(np.float32)
        transltransf = mat[:, 2].astype(np.float32)
        return affineTransform(lintransf, transltransf)
    
    def accumulate(self, affinetransform2):
        lintransf_res = np.matmul( self.lintransf, affinetransform2.lintransf)
        transltransf_res = self.transltransf + affinetransform2.transltransf
        return affineTransform(lintransf_res, transltransf_res)
    
    def getAffTransfAs3x2mat(self):
        mat = np.column_stack((self.lintransf,self.transltransf))
        return mat
    
    def addTransl(self, transl):
        self.transltransf+=transl

    def applyToPoint(self,point):
        #point being (x,y)
        x0,y0= point

        x1 = self.lintransf[0,0] * x0 + self.lintransf[0,1] * y0 + self.transltransf[0]
        y1 = self.lintransf[1,0] * x0 + self.lintransf[1,1] * y0 + self.transltransf[1]

        return x1,y1
    
    def copy(self):
        return affineTransform(self.lintransf, self.transltransf)


def align_stack(image_sequence, method, callbk_tick_fn=None):
    '''
    Aligns a 3D stack along the z-axis using the SIFT aligment method.

    Parameters:
        image_sequence: a 3D numpy array with the images (greyscale) to be aligned along the X plane

    method: a dictionary that sets the alignement method, like following:
        'translation':True,
        'rotation':False, 'shearing_x':True, 'shearing_y':False,
        'scaling':False, 'stretching_x':False, 'stretching_y':True, 
        'affine':False 

        'translation' Can be overriden by affine setting below
        Either 'rotation', 'shearing_x' or 'shearing_y' . If rotation is true the shearings are ignored
        Either 'scaling', 'stretching_x' or 'stretching_y'. If scaling is true, the stretchings are ignored
        'affine'. This setting mean that all parameters in the affine matrix are free to be optimized for alignement.
                If affine is true all the others are ignored.
    
    callbk_tick_fn: call back function that will be called between iterations

    Returns:
        The 3D stack aligned data.
        It is likely to be larger in pixels along the XY plane to compensate for the translation and other transforms

    Estimated number of callbacks = (nslices-1) + 3 + nslices = 2*nslices +2
        
    '''

    #prepare list of transformation matrice and canvas 
    aff_transf_0 = affineTransform.from2x3mat(np.float32([[1, 0, 0], [0, 1, 0]])) #no transform, identity and no translation
    mats = [aff_transf_0]

    #match each image with the previous image as reference
    for i in range(len(image_sequence) - 1):
        ref = image_sequence[i]
        new = image_sequence[i + 1]
        print(f'Matching slice {i + 1} and slice {i + 2} ...')    
        
        matching_failed=False
        #feature matching
        matches = featureMatchingDS(ref, new)
        #print(f"len(matches): {len(matches)}")
        #print(f"len(matches[0]): {len(matches[0])}")
        if len(matches[0])>0:
            x_ref = matches[0][:, 0, 0]
            y_ref = matches[0][:, 0, 1]
            coor_ref = [x_ref, y_ref]

            x_new = matches[1][:, 0, 0]
            y_new = matches[1][:, 0, 1]
            coor_new = [x_new, y_new]
        
            #find optimal transformation
            t0, bounds = findInitialandBounds(method)
            res = least_squares(residual(coor_ref, coor_new, method), t0, bounds=bounds)
            print('--------------------------------------------------------------------')
            print('Optimisation result: ')
            print('Cost: ', res.cost)
            print('Message: ', res.message)
            print('Success: ', res.success)
            print('Parameters: ', res.x)
        
            if not res.success:
                matching_failed=True
        else:
            matching_failed=True

        if not matching_failed:
            T, affine, translate = findMatrix(method, res.x)
            T1 = affineTransform.from2x3mat(T)
        else:
            T1= aff_transf_0.copy() #Use the identity
            #Use last one calculated
            #T1 = mats[-1].copy()

        mats.append(T1)

        if not callbk_tick_fn is None:
            callbk_tick_fn()
        
    print('** Completed getting the slice-to-slice transforms.** ')

    print('Accumulating transforms along the stack')
    mat_accum = [aff_transf_0] #First image transform
    mat0 = aff_transf_0
    for i in range(1,len(mats)):
        mat1 = mat_accum[i-1].accumulate(mats[i]) #Calculate new matrix
        mat_accum.append(mat1)
        mat0=mat1
    
    print('Completed calculation of accumalation transform matrices.')

    if not callbk_tick_fn is None:
            callbk_tick_fn()

    canvas = []
    #output = []

    for T in mat_accum:
        #calculate image size so that all images are completely on the canvas
        width_new, length_new = new.shape
        width_ref, length_ref = ref.shape
        width = max(width_new, width_ref)
        length = max(length_new, length_ref)
        
        before_points = np.float32([[0, 0], [0, width], [length, width], [length, 0]]).reshape(-1, 1, 2)
        #after_points = cv.transform(before_points, T)
        after_points = cv.transform(before_points, T.getAffTransfAs3x2mat() )
        
        list_of_points = np.float32([[0, 0], [0, width], [length, width], [length, 0]]).reshape(-1, 1, 2)
        list_of_points = np.append(list_of_points, after_points, axis=0)
        
        [x_minT, y_minT] = np.int32(list_of_points.min(axis=0).ravel() - 0.5)
        [x_maxT, y_maxT] = np.int32(list_of_points.max(axis=0).ravel() + 0.5)
        
        canvas.append([x_minT, y_minT, x_maxT, y_maxT])
        
    x_min = min(canvas[i][0] for i in range(len(canvas)))
    y_min = min(canvas[i][1] for i in range(len(canvas)))
    x_max = max(canvas[i][2] for i in range(len(canvas)))    
    y_max = max(canvas[i][3] for i in range(len(canvas)))
    
    if not callbk_tick_fn is None:
            callbk_tick_fn()

    ysize = y_max - y_min
    xsize = x_max - x_min

    translation_distT = [-x_min, -y_min]
    print(f'translation_distT:{translation_distT} , need to be added to the transform')
    

    
    for i0 in range(len(mat_accum)):
        mat_accum[i0].addTransl(translation_distT)

    if not callbk_tick_fn is None:
        callbk_tick_fn()

    print('Beggining applying transformations to each image in the stack')

    image_sequence_res = np.zeros( (len(image_sequence), ysize,xsize), dtype=image_sequence.dtype)

    for i in range(len(image_sequence)):
        mat0 = mat_accum[i].getAffTransfAs3x2mat()
        image_t = cv.warpAffine(image_sequence[i], mat0, (xsize, ysize))
        #plt.imsave(f'../sequence/{i + 1}.jpg', image_t, cmap='gray')
        #output.append(image_t)
        image_sequence_res[i,:,:] = image_t

        if not callbk_tick_fn is None:
            callbk_tick_fn()

    print('Alignment completed.')

    #TODO: Consider returning translations and alignment matrices used between slices in
    # case user wants to plot

    return image_sequence_res


