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


import numpy as np

def fft_bandpass_filter_dirx_2D(data2d , freq_low= 0, freq_high = 5):
    '''
    A simple directional bandpass filter
    Parameters:
        data2d : two-dimensional numpy array with image to be filtered
        freq_low, freq_high: lower and upper bounds of the banbpass filter in 1/pixel units

    '''
    # Applies a FFT filtering along the x-axis (axis=1)
    # It cuts any signal with frequencies along the x-axis, with values between freq_low and freq_high
    dataret = None
    print(f"data2d.shape: {data2d.shape}")
    
    if data2d.ndim ==2:
        data_dir_fft = np.fft.rfft(data2d, axis=1)
        #print(f"data_dir_fft.shape: {data_dir_fft.shape}")

        data_dir_fft[:,freq_low:freq_high] = 0

        dataret = np.fft.irfft(data_dir_fft, axis=1)
    
    return dataret

