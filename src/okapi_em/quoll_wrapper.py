import numpy as np
import quoll.frc
import tempfile
import tifffile

#Quoll (FRC) as an optional package?
is_quoll_available = True
try:
    from quoll.frc import oneimg
    from quoll.io import reader, tiles
except:
    is_quoll_available=False



def getFRC(data2square):
    '''
    Calculates the FRC of the whole image

    '''
    res=None

    # TODO: check data is square. If not, we may need to crop centre

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

def getTiledFRC(data2d):
    print("Not implemented yet")
    pass
