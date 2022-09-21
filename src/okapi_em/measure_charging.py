import math

#import matplotlib.pyplot as plt
import napari
import numpy as np
from skimage import io


def open_image(input_image):
    """Opens an image as np array

    Parameters
    ----------
    input_image : str, path-like
        Path to input image

    Returns
    -------
    np array
        Numpy array containing the image.
    """
    img = io.imread(input_image)
    return np.array(img)


def extract_tiles(img, tile_size):
    """Split image `img` into square tiles of equal size `tile_size`

    The image is zero-padded to ensure that the entire image is sampled.

    Parameters
    ----------
    img : np array
        Numpy array containing the image
    tile_size : int
        Length of each side of the square tile

    Returns
    -------
    list
        List containing numpy arrays of shape `(tile_size, tile_size)`
        taken from the image. 
    """
    # How many tiles will we use
    nTilesX = math.ceil(img.shape[1] / tile_size)
    nTilesY = math.ceil(img.shape[0] / tile_size)

    # Calculate padded size and create empty image
    paddedX = math.ceil(nTilesX * tile_size)
    paddedY = math.ceil(nTilesY * tile_size)
    padded_img = np.zeros(shape=(paddedY, paddedX))

    # Place image in the centre of the padded image
    ybounds = [int((paddedX/2) - (img.shape[1]/2)), int((paddedX/2) + (img.shape[1]/2))]
    xbounds = [int((paddedY/2) - (img.shape[0]/2)), int((paddedY/2) + (img.shape[0]/2))]
    padded_img[xbounds[0]:xbounds[1], ybounds[0]:ybounds[1]] = img

    # Split into tiles
    tiles = []
    for i in range(nTilesX):
        for j in range(nTilesY):
            tile_xbounds = [j*tile_size, (j*tile_size) + tile_size]
            tile_ybounds = [i*tile_size, (i*tile_size) + tile_size]
            tile = padded_img[tile_xbounds[0]:tile_xbounds[1], tile_ybounds[0]:tile_ybounds[1]]
            tiles.append(tile)
    return tiles, (nTilesX, nTilesY)


def calc_2D_fft(img):
    """Calculate and return the magnitude of the 2D-FFT

    First the 2D-Fourier transform is calculated, and the zero-frequency
    component is shifted to the centre of the spectrum. The absolute value
    of the shifted FT is returned.

    Parameters
    ----------
    img : np array
        Numpy array containing the image

    Returns
    -------
    np array
        2D Numpy array containing the magnitude of the 2D Fourier transform.
    """
    ft = np.fft.fft2(img)
    ft_shifted = np.fft.fftshift(ft)
    ft_magnitudes = 20*np.log(np.abs(ft_shifted))
    return ft_magnitudes


def calc_Q(ft, horizontal=True):
    """Calculate charging artefacts factor, `Q`.
    
    High `Q` means more charging artefacts are occuring in the patch.

    Q = max of the column-averaged FT / median of the column-averaged FT

    If charging artefacts are vertical, then the row-averaged FT is used.

    Parameters
    ----------
    ft : np array
        2D Numpy array containing the magnitude of the 2D FT
    horizontal : bool
        Direction of charging artefacts. Default value is True.
        If artefacts are vertical, set to False.

    Returns
    -------
    Q
        Charging artefacts factor
    """
    if horizontal == True:
        average = np.average(ft, axis=0)
    else:
        average = np.average(ft, axis=1)

    Q = np.max(average) / np.median(average)
    return Q


def reassemble_tiles(tiles, nTiles):
    """Reassembles tiles into whole images

    Parameters
    ----------
    tiles : list
        List containing 2D numpy arrays of the tiles.
    nTiles : tuple
        (Number of tiles in x, number of tiles in y)

    Returns
    -------
    np array
        2-D np array containing the whole image put together
    """
    nTilesX, nTilesY = nTiles
    tile_size = tiles[0].shape[0]
    img = np.zeros(shape=(nTilesY * tile_size, nTilesX * tile_size))
    counter = 0
    for i in range(nTilesX):
        for j in range(nTilesY):
            tile_xbounds = [j*tile_size, (j*tile_size) + tile_size]
            tile_ybounds = [i*tile_size, (i*tile_size) + tile_size]
            img[tile_xbounds[0]:tile_xbounds[1], tile_ybounds[0]:tile_ybounds[1]] = tiles[counter]
            counter += 1
    return img


def generate_heatmap(img, tile_size, horizontal=True):
    """Generates a heatmap of the charging artefacts

    Parameters
    ----------
    img : np array
        2D np array of the image
    tile_size : int
        Size of tiles to use to calculate Q
    horizontal : bool, optional
        direction of charging artefacts, by default True

    Returns
    -------
    np array
        Heatmap that can be overlaid on the original image
        showing extent of charging artefacts
    """
    tiles, (nTilesX, nTilesY) = extract_tiles(img, tile_size)
    Qs = [calc_Q(calc_2D_fft(tile)) for tile in tiles]
    Q_overlay = [np.full(shape=tiles[0].shape, fill_value=Q) for Q in Qs]
    overlay = reassemble_tiles(Q_overlay, (nTilesX, nTilesY))
    return overlay    


# if __name__ == "__main__":
#     input_image = "2022-02-22_20_12_44_1_18_mode_0.tif"
#     img = open_image(input_image)
#     tiles, (nTilesX, nTilesY) = extract_tiles(img, 512)
#     overlay = generate_heatmap(img, 256)

#     # plt.imshow(img2, cmap="Greys_r")
#     # plt.imshow(overlay, cmap="viridis", alpha=0.3)
#     # plt.colorbar()
#     # plt.show()

#     nv = napari.Viewer()
#     nv.add_image(img, name="Original image")
#     nv.add_image(overlay, name="Charging artefacts", opacity=0.3, colormap="viridis")
#     napari.run()
#     print("Done")
