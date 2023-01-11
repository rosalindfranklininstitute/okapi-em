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

'''
A wrapper to help using quoll for FRC calculations
'''

import numpy as np
from inspect import getmembers, isfunction

# Quoll (FRC) as an optional package?
is_quoll_available = True
try:
    from quoll.frc import oneimg
    from quoll.io import reader, tiles
    from quoll.frc import frc_calibration_functions as cfuncs
except:
    is_quoll_available=False

def getTiledFRC(data2d, cf_in, tilesize=256, pixel_size_nm=1.0):
    '''
    Calculates the FRC of the image by splitting it into square tiles of shape 
    (tilesize,tilesize)

    Params:
    data2d: np.array containing image data
    tilesize: int, min 128. Length of one side of the square tile in pixels.
    pixel_size_nm: float, >0. Physical size of each pixel.
    cf: str. Name of the function from frc_calibration_functions to use.

    Returns: a tuple with two elements
        dfres: a list of FRC values calculated of each tile (no coordinates though)
        FRC heatmap: An image with same dimensions of the original image but as a tiled heatmap of the
            FRC value in each of the tiles

    '''

    res = None

    #check data is square. If not, we may need to crop centre

    if data2d.ndim == 2:
        
        # load image as a Quoll Image object
        Img = reader.Image(
            img_data=data2d,
            pixel_size=pixel_size_nm,
            unit="nm"
        )

        # Get calibration function
        if cf_in == "None":
            cf = None
            print("No calibration performed")
        else:
            cf_list = dict(getmembers(cfuncs, isfunction))
            cf = cf_list[cf_in]
            print(f"Calibration function is {cf_in}: {cf}")

        # Calculate the FRC
        dfres = oneimg.calc_local_frc(
            Img,
            tile_size=tilesize,
            calibration_func=cf
        )

        res = dfres, get_resolution_heatmap(Img, dfres)
    
    else:
        raise ValueError("Data is not 2D.")

    return res

def get_resolution_heatmap(
    Image,
    results_df
):
    """Generate a heatmap of resolutions
    
    Optionally display the heatmap overlaid on the original image,
    save this overlay, or save the heatmap alone.

    Args:
        Image (reader.Image): Quoll.reader.Image object
        results_df (pd.DataFrame): output of `oneimg.calc_local_frc`

    Returns:
        _type_: _description_
    """

    tileshape = list(Image.tiles.values())[0].shape
    restiles = [np.full(shape=tileshape, fill_value=res) for res in np.array(results_df.Resolution)]
    heatmap =  tiles.reassemble_tiles(
        tiles=restiles,
        nTiles=Image.tile_arrangement,
    )
        
    return heatmap


