def fourier_transform(X,Y):
    '''
    '''

    T = X[1] - X[0]
    N = len(X)

    yf = fftpack.fft(Y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    amplitudes = 2.0/N * np.abs(yf[0:N//2])

    return xf, amplitudes

def fourier_amplitude(X,Y, cutoffs = [0.0005,0.1]):
    '''
    Get amplitudes from X,Y data via fourier transform.

    Parameters
    ----------
    X: 1D numpy array
        array of x values.
    Y: 1D numpy array
        array of y values.
    cutoff: list of floats
        lower and upper frequencies (in Hz) in which to narrow the area to be
        searched for maximum amplitudes.

    Returns
    -------
    np.nan_to_num(max_amp): 1D numpy array
        An array of maximum amplitude values within the set cutoff frequency
        bounds.
    '''

    T = X[1] - X[0]
    N = len(X)

    yf = fftpack.fft(Y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    amplitudes = 2.0/N * np.abs(yf[0:N//2])

    inds = np.where((xf>cutoffs[0])&(xf<cutoffs[1]))[0]

    amplitudes = amplitudes[inds]

    max_amp = np.amax(amplitudes)

    return np.nan_to_num(max_amp)

def Peak_finder(intensity, thres=0.1, min_dist=1, max_inten = 1e100, min_inten = -1e100, d_type = 'null'):

    '''Peak detection routine.
    Modified from PeakUtils:
    https://github.com/atjacobs/PeakUtils/tree/master/peakutils
    Finds the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks.

    Parameters
    ----------
    y : ndarray
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    max_inten: float
        peaks will not be detected above this threshold.
    min_inten: float
        peaks will not be detected below this threshold.
    d_type: str
        not used

    Returns
    -------
    ndarray
        Array containing the indexes of the peaks that were detected
    '''

    test = np.where((intensity < max_inten))[0]
    thres *= np.max(intensity[test]) - np.min(intensity[test])

    # find the peaks by using the first order difference
    diff = np.diff(intensity)

    peaks_indices = np.where((np.hstack([diff, 0.]) < 0.)
                     & (np.hstack([0., diff]) > 0.)
                     & (intensity > thres) & (intensity < max_inten) & (intensity > min_inten))[0]

    if peaks_indices.size > 1 and min_dist > 1:
        highest = peaks_indices[np.argsort(intensity[peaks_indices])][::-1]
        rem = np.ones(intensity.size, dtype=bool)
        rem[peaks_indices] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peaks_indices] = False

        peaks_indices = np.arange(intensity.size)[~rem]

    # find beginning of peaks
    peak_starts = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]-2
        while diff[cursor] > 0:
            cursor -= 1
        peak_starts.append(cursor)

    # find end of peaks
    peak_ends = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]+2
        if int(cursor) >= len(diff):
            cursor = len(diff)-1

        while diff[cursor] < 0:
            cursor += 1
            if int(cursor) >= len(diff):
                cursor = len(diff)-1
                break
        peak_ends.append(cursor)

    return {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}
