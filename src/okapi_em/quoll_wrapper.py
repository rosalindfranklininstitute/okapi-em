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

        results_df = oneimg.calc_frc_res(Img)

        tempdir.cleanup()

        print (f"results_df : {results_df}")
        res= results_df


    return res

def getTiledFRC(data2d):
    print("Not implemented yet")
    pass
