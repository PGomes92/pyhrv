def spectrogram(nni=None, rpeaks=None, nfft=2**16):
    nn = tools.check_input(nni, rpeaks)

    # Resampling (with 4Hz) and interpolate
    # Because RRi are unevenly spaced we must interpolate it for accurate PSD estimation.
    fs = 4
    t = np.cumsum(nn)
    t -= t[0]
    f_interpol = sp.interpolate.interp1d(t, nn, 'cubic')
    t_interpol = np.arange(t[0], t[-1], 1000. / fs)
    nn_interpol = f_interpol(t_interpol)

    # Adapt 'nperseg' according to the total duration of the NNI series (5min threshold = 300000ms)
    if t.max() < 300000:
        nperseg = nfft
    else:
        nperseg = 300

    f, t, Sxx = sp.signal.spectrogram(
        x=nn_interpol,
        fs=fs,
        window='hamming',
        nperseg=nperseg,
        scaling='spectrum'
    )

    # Plot of Spectrogram
    fig, ax = plt.subplots(figsize=(25, 5))
    ax.pcolormesh(t, f, Sxx)
    ax.plot([t[0], t[-1]], [0, 0.3])
    ax.set_ylabel('Frequency (Hz)')
    ax.legend()
    ax.set_ylim([0, 0.4])

    # X-Axis configuration
    # Set x-axis format to seconds if the duration of the signal <= 60s
    if t[-1] <= 60:
        ax.set_xlabel('Time [s]')
    # Set x-axis format to MM:SS if the duration of the signal > 60s and <= 1h
    elif 60 < t[-1] <= 3600:
        ax.set_xlabel('Time [MM:SS]')
        formatter = mpl.ticker.FuncFormatter(lambda ms, x: str(dt.timedelta(seconds=ms))[2:])
        ax.xaxis.set_major_formatter(formatter)
    # Set x-axis format to HH:MM:SS if the duration of the signal > 1h
    else:
        ax.set_xlabel('Time [HH:MM:SS]')
        formatter = mpl.ticker.FuncFormatter(lambda ms, x: str(dt.timedelta(seconds=ms)))
        ax.xaxis.set_major_formatter(formatter)

    plt.show()