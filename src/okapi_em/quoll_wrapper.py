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
import tempfile
import tifffile

#Quoll (FRC) as an optional package?
is_quoll_available = True
try:
    #import quoll.frc
    from quoll.frc import oneimg
    from quoll.io import reader, tiles
except:
    is_quoll_available=False



def getFRC(data2square):
    '''
    Calculates the FRC of the whole image

    Note that this uses the miplib function directly, which has the mibplib's author calibration.
    Avoid using this function. This is just sitting here as it may be usefule in the future.

    '''
    res=None

    # check data is square. If not, we may need to crop centre

    if data2square.ndim ==2:

        #Create temporary file, because quoll+miplib only works with quoll.io.reader.Img objects,
        # and can only be created from files
        
        tempdir = tempfile.TemporaryDirectory()
        
        image_filename = tempdir.name+'/image_temp.tif'
        print(f'temporary image_filename: {image_filename}')

        tifffile.imsave(image_filename,data2square)

        Img = reader.Image(
            filename=image_filename
        )
        Img.pixel_size=1

        #results_df = oneimg.calc_frc_res(Img) #Images must be square
        try:
            import miplib.ui.cli.miplib_entry_point_options as opts
            from miplib.data.io import read as miplibread

            miplibImg = miplibread.get_image(Img.filename)
            miplibargs = opts.get_frc_script_options([None])

            results_frc1 = oneimg.miplib_oneimg_FRC_calibrated(miplibImg, miplibargs)
            print (f"results_df : {results_frc1}")


            frc_value = results_frc1.resolution["resolution"]
            res= frc_value
        except:
            print("Error using miplib")

        tempdir.cleanup()
        
    return res


def getTiledFRC(data2d, tilesize=256):
    '''
    Calculates the FRC of the image by splitting it into square tiles of shape (tilesize,tilesize)

    Returns: a tuple with two elements
        dfres: a list of FRC values calculated of each tile (no coordinates though)
        FRC heatmap: An image with same dimensions of the original image but as a tiled heatmap of the
            FRC value in each of the tiles 

    '''

    res=None

    #check data is square. If not, we may need to crop centre

    if data2d.ndim ==2:

        tempdir = tempfile.TemporaryDirectory()
            
        image_filename = tempdir.name+'/image_temp.tif'
        print(f'temporary image_filename: {image_filename}')

        tifffile.imsave(image_filename,data2d)

        Img = reader.Image(
                filename=image_filename
            )
        Img.pixel_size=1 #prevent errors. quoll uses this to estimate resolution in physical units

        dfres = oneimg.calc_local_frc(Img, tilesize, tempdir.name)

        res = dfres,get_resolution_heatmap(Img, dfres)

        tempdir.cleanup()

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


