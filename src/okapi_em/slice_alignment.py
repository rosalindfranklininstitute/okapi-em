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
    'translation':True,
    'rotation':False, 'shearing_x':True, 'shearing_y':False, # If rotation is true the shearings are ignored
    'scaling':False, 'stretching_x':False, 'stretching_y':True, # If scaling is true, the stretchings are ignored
    'linear':False # If linear is true all the others are ignored. This setting mean that all parameters in the linear matrix are free to adjust in the optimization
}


def featureMatchingDS(img1, img2):
    # Initiate SIFT detector
    sift = cv.SIFT_create()
    # find the keypoints and descriptors with SIFT
    kp1, des1 = sift.detectAndCompute(img1,None)
    kp2, des2 = sift.detectAndCompute(img2,None)
    #This is probably unnecessary to do it like this
    # In a stack, most images will have SIFT computed twice
    # TODO: Create a class that is initialised with image stack
    # and precalculates Sift's keypoints and descriptors

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
    #img_match_gauss = cv.drawMatchesKnn(img1,kp1,img2,kp2,matched_filtered,None,flags=cv.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS)
    #plt.figure(figsize = (20,20))
    #plt.imshow(img_match_gauss)
    
    return reference_pts, new_pts


def findMatrix(method, params):

    #Default to no transform, no translation
    T= np.array( [[1,0, 0],
                  [0,1, 0]], dtype='float32')
    #T is the transform matrix
    #T = [a0 a1 | t0]
    #    [a2 a3 | t1]

    #print(f"findMatrix, params:{params}")

    if method['translation']==True:
        T += np.array([ [0,0,params[0]],
                        [0,0,params[1]] ], dtype='float32')

    if method['linear'] == True:
        T += np.array([[params[2], params[3], 0],
                       [params[4], params[5], 0]], dtype='float32')
    else:
        if method['rotation'] == True:
            T += np.array([[np.cos(params[2]), - np.sin(params[2]), 0], 
                           [np.sin(params[2]), np.cos(params[2]), 0]], dtype='float32')
        elif method['shearing_x']==True:
            T += np.array([[ 0 , params[2] , 0], 
                           [0 , 0 , 0]], dtype='float32')
        elif method['shearing_y']==True:
            T += np.array([[ 0 , 0 , 0], 
                           [params[2] , 0 , 0]], dtype='float32')

        if method['scaling']==True:
            T *= np.array([[ params[3] , params[3] , 1], 
                           [params[3] , params[3] , 1]], dtype='float32')
        elif method['stretching_x']==True:
            T *= np.array([[ params[3] , 1 , 1], 
                           [params[3] , 1 , 1]], dtype='float32')
        elif method['stretching_y']==True:
            T *= np.array([[ 1 , params[3] , 1], 
                           [1 , params[3] , 1]], dtype='float32')

    if T is None:
        raise ValueError("No suitable method transform found")
            
    linearT = T[:, :2]
    translateT = T[:, 2]
    return T, linearT, translateT


def transform(method, param, points):
    """
    Applies a 2D affine transformation (linear+translate) to points (x,y)

    This function uses findMatrix method to get the transformation before applying to points
    
    Params:
        method: a dictionary with boolean entries for
        translation, linear, rotation, shearing_x, shearing_y,
        scaling, stretching_x, stretching_y

        param: 1D array with elements being the parameters of the transformation
        indexes 0 and 1 are x and y translation
        Rest of the indexes depend on the transformation selected with 'method'
    """
    
    # Use findMatrix code
    T, linearT,translateT = findMatrix(method, param)

    myAff = cAffineTransform(linearT,translateT)
    #Apply matrix to points
    p_res=myAff.applyToPoint(points)
    return p_res


def residual(ref, new, method):
    def r(t):
        rx = transform(method, t, new)[0] - ref[0]
        ry = transform(method, t, new)[1] - ref[1]
        return np.append(rx, ry)
    return r

def setInitialandBounds(method, margin_t=200, margin_a=0.01):
    """ Sets initial vallues of params depending on the optimization required and bounds
    margin_t sets that variation allowed in translation
    margin_a sets the variation allowed with other parameters around 0 or around 1, depending on the parameter

    Returns:
        params_init: initial parameters as a 1D numpy array
        bounds: respective bounds as an array with 2 rows, with 1st row being
            the lower bounds and second row being the upper bounds.
    """
    
    params_init=np.array([0, 0,   1, 0, 0, 1]) #tr_x, tr_y, rot/shear/p2 , scale/stretch/p3, p4, p5
    bounds = np.array([ [-0.0001, -0.0001, 0.9999, -0.0001, -0.0001, 0.9999],
                       [ 0.0001,  0.0001, 1.0001,  0.0001,  0.0001, 1.0001] ])
    
    if method['translation']==True:
        #First two parameters are the translation
        bounds[:,:2] = np.array([ [-margin_t, -margin_t],
                                 [ margin_t,  margin_t] ])
    
    if method['linear']==True:
        #Linear part of the array is free to vary
        bounds[:,2:] = np.array([ [1-margin_a, -margin_a, -margin_a, 1-margin_a],
                                 [1+margin_a,  margin_a,  margin_a, 1+margin_a]])
        
    else:
        #No p4 or p5 needed
        params_init=params_init[:-2]
        bounds=bounds[:,:-2]

        if method['scaling']==False and method['stretching_x']==False and method['stretching_y']==False:
            #reduce the parameters and bounds
            # no p3 needed
            params_init=params_init[:-1]
            bounds=bounds[:,:-1]

            if method['rotation']==False and method['shearing_x']==False and method['shearing_y']==False:
                #no p2 needed
                params_init=params_init[:-1]
                bounds=bounds[:,:-1]
        else:
            if method['scaling']==True or method['stretching_x']==True or method['stretching_y']==True:
                params_init[3] = 1 #variations of scalings around one
                bounds[0,3] = 1-margin_a
                bounds[1,3] = 1+margin_a
            
            if method['rotation']==True or method['shearing_x']==True or method['shearing_y']==True:
                params_init[2] = 0 #variations of angle and shearing around zero
                bounds[0,2] = -margin_a
                bounds[1,2] = +margin_a
    
    #print(f"params_init:{params_init}")
    #print(f"bounds:{bounds}")

    return params_init, bounds

class cAffineTransform():
    '''
    Class for affine transforms,
    typically a matrix for rotation and/or skewing
    and a translation.

    Contains functions to do operations with these matrices
    '''
    def __init__(self,lintransf, transltransf):
        '''
        Args:
            lintransf: a 2x2 numpy matrix with the rotation/skewing values
            transltransf: translation values in matrix with 2 in format [x,y] 
        
        '''
        self.lintransf = lintransf.copy()
        self.transltransf= transltransf.copy()
    
    @staticmethod
    def from2x3mat(mat):
        '''
        Creates a affineTransform object from a 2x3 matrix. It is a constructor type of function.
        Args:
            mat: a 2x3 matrix with the linear transform in the first two columns
                and the translation in the last column
        Returns:
            a new affineTransform object
        '''
        lintransf = mat[:, :2].astype(np.float32)
        transltransf = mat[:, 2].astype(np.float32)
        return cAffineTransform(lintransf, transltransf)

    def getAffTransfAs3x2mat(self):
        '''
        Gets a 2x3 matrix numpy of this affineTranform object.
        Returns:
            a numpy array 2x3, with the linear transform in the first two columns
                and the translation in the last column
        '''
        mat = np.column_stack((self.lintransf,self.transltransf))
        return mat

    def accumulate(self, affinetransform2):
        '''
        Accumulates affine transforms.
        This is not a multiplication or addition.
        In practice, the linear transform is multiplied and the translation is added.
        Args:
            affinetransform2, affine transform to accumulate
        '''
        lintransf_res = np.matmul( self.lintransf, affinetransform2.lintransf)
        transltransf_res = self.transltransf + affinetransform2.transltransf
        return cAffineTransform(lintransf_res, transltransf_res)
    

    def addTransl(self, transl):
        '''
        Adds a translation [x,y] to the this affineTransform object
        '''
        self.transltransf+=transl

    def applyToPoint(self,point):
        '''
        Applies this affineTransformation to a point/vector
        Args:
            point: coordinates in format (x,y)

        Returns:
            a new point (xr,yr) that is the tranformed of [x,y]
        '''
        #point being (x,y)
        x0,y0= point

        x1 = self.lintransf[0,0] * x0 + self.lintransf[0,1] * y0 + self.transltransf[0]
        y1 = self.lintransf[1,0] * x0 + self.lintransf[1,1] * y0 + self.transltransf[1]

        return x1,y1
    
    def copy(self):
        '''
        Makes a copy of itself
        '''
        return cAffineTransform(self.lintransf, self.transltransf)


def align_stack(image_sequence, method, callbk_tick_fn=None):
    '''
    Aligns a 3D stack along the z-axis using the SIFT aligment method.

    Args:
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
    aff_transf_0 = cAffineTransform.from2x3mat(np.float32([[1, 0, 0], [0, 1, 0]])) #no transform, identity and no translation
    mats = [aff_transf_0]

    #match each image with the previous image as reference
    for i in range(len(image_sequence) - 1):
        ref = image_sequence[i]
        new = image_sequence[i + 1]
        print("---------------------------------------")
        print(f'Matching slice {i} and slice {i+1} ...')    
        
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
            t0, bounds = setInitialandBounds(method)
            res = least_squares(residual(coor_ref, coor_new, method), t0, bounds=bounds)
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
            T1 = cAffineTransform.from2x3mat(T)
        else:
            T1= aff_transf_0.copy() #Use the identity
            #Use last one calculated
            #T1 = mats[-1].copy()

        mats.append(T1) #Stores affine transform, it will be applied later

        if not callbk_tick_fn is None:
            callbk_tick_fn()

    print("Completed")
    print("Accumulating transforms along the stack")
    mat_accum = [aff_transf_0] #First image transform
    mat0 = aff_transf_0
    for i in range(1,len(mats)):
        mat1 = mat_accum[i-1].accumulate(mats[i]) #Calculate new matrix
        mat_accum.append(mat1)
        mat0=mat1
    
    print('Completed calculation of accumulation transform matrices.')

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
        after_points = cv.transform(before_points, T.getAffTransfAs3x2mat() ) #  we may not need the cv.transform but we can maybe use the applyToPoint function
        
        #LMAP: I am not really sure what is going on below
        # It appears that list_of_points is being assigned values like before_points
        # and then a new dimension added to include after_points
        # Nevertheles: it appears to be working, for the next step
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


class cFeatureMatching():
    """
    Class to handle feature matching.
    It tries to be more efficient than the function based method by precalculating
    the SIFT key points and descriptors for the volumes
    """

    def __init__(self, image_sequence):
        #TODO: Check image is valid and ok to run with CV2, otherwise convert

        self.image_sequence = image_sequence

        self._sift = cv.SIFT_create()

        #Precalculate kp and ds for all image sequences
        self.kps = []
        self.dss = []

        for im0 in image_sequence():
            kp0, ds0 = self._sift.detectAndCompute(im0,None)
            self.kps.append(kp0)
            self.dss.append(ds0)
    
    def featureMatchingDS(self, idx0, idx1, d_to_dmedian_perc_filt = 0.1, rel_dist_to_second_match=0.65):
        """
        idx0, idx1 are indexes of the images to compare
        Returns matching points in first and second image
        using criteria set with d_to_dmedian_perc_filt and rel_dist_to_second_match

        TODO: explain parameters, d_to_dmedian_perc_filt and rel_dist_to_second_match
        """
        
        #Check indices are valid
        if idx0>=len(self.kps1) or idx0<0 or idx1>len(self.kps1) or idx1<0 or idx0==idx1:
            raise ValueError(f"idx0:{idx0} or idx1:{idx1} are not valid for number of entries:{len(self.kps1)}.")
        
        if rel_dist_to_second_match<=0 or rel_dist_to_second_match>=1:
            raise ValueError("rel_dist_to_second_match is not valid")
        
        #Do matching using opencv BFMatcher
        bf = cv.BFMatcher() #Brute Force description matcher using default settings
        
        ds0= self.dss[idx0]
        ds1 = self.dss[idx1]
        matches = bf.knnMatch(ds0,ds1,k=2) #Gets 2 matches for each query descriptor

        #Filters by descriptor distance
        # Filters best matches by distance by its certainty by comparing with second best match distance
        # If there is another match that is similar, don't take it into account
        good_no_list = []
        for m, n in matches:
            if m.distance < rel_dist_to_second_match*n.distance:
                #good.append([m])
                good_no_list.append(m)
        
        kp0 = self.kps[idx0]
        kp1 = self.kps[idx1]

        #Gets coordinates of filtered points in first and second image
        coords0 = np.float32([kp0[m.queryIdx].pt for m in good_no_list])
        coords1 = np.float32([kp1[m.trainIdx].pt for m in good_no_list])

        diff_vects = coords1-coords0

        slope_diff = np.arctan2(diff_vects[:,1], diff_vects[:,0])

        dist_diff = np.linalg.norm(diff_vects, axis=1)
        dist_diff_median = np.median(dist_diff)

        dist_select = np.logical_and(dist_diff > dist_diff_median * (1.0-d_to_dmedian_perc_filt), 
                                     dist_diff < dist_diff_median * (1.0+d_to_dmedian_perc_filt))
        
        #Filter by the intensity of the points
        coords0_int = coords0.astype(int)
        img_values0 = [self.image_sequence[idx0][x, y] for y, x in coords0_int]
        coords1_int = coords1.astype(int)
        img_values1 = [self.image_sequence[idx1][x, y] for y, x in coords1_int]

        #Create mask by looking at pixel values between points need to have intensity >1
        points_select = np.logical_and(np.array(img_values0) > 1, np.array(img_values1) > 1)

        #Selection mask based from the two filters dist_select and points_select
        select = dist_select * points_select

        # matched_filtered = np.array(good_no_list)[select] #not needed unless plotting
        coords0_filt = coords0[select]    
        coords1_filt = coords1[select]

        print(f'{len(coords0_filt)} good matches found')

        return coords0_filt, coords1_filt