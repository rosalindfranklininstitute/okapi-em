import numpy as np

def fft_bandpass_filter_dirx_2D(data2d , freq_low= 0, freq_high = 5):
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

