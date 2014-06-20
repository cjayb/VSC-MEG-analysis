import numpy as np
from mne.viz import _prepare_topo_plot, plot_topomap, _check_delayed_ssp, _draw_proj_checkbox

COLORS = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#473C8B', '#458B74',
          '#CD7F32', '#FF4040', '#ADFF2F', '#8E2323', '#FF1493']

DEFAULTS = dict(color=dict(mag='darkblue', grad='b', eeg='k', eog='k', ecg='r',
                           emg='k', ref_meg='steelblue', misc='k', stim='k',
                           resp='k', chpi='k', exci='k', ias='k', syst='k'),
                units=dict(eeg='uV', grad='fT/cm', mag='fT', misc='AU'),
                scalings=dict(eeg=1e6, grad=1e13, mag=1e15, misc=1.0),
                scalings_plot_raw=dict(mag=1e-12, grad=4e-11, eeg=20e-6,
                                       eog=150e-6, ecg=5e-4, emg=1e-3,
                                       ref_meg=1e-12, misc=1e-3,
                                       stim=1, resp=1, chpi=1e-4, exci=1,
                                       ias=1, syst=1),
                ylim=dict(mag=(-600., 600.), grad=(-200., 200.),
                          eeg=(-200., 200.), misc=(-5., 5.)),
                titles=dict(eeg='EEG', grad='Gradiometers',
                            mag='Magnetometers', misc='misc'))

def plot_evoked_topomap(evoked, times=None, ch_type='mag', layout=None,
                        vmax=None, vmin=None, cmap='RdBu_r', sensors='k,',
                        colorbar=True, scale=None, scale_time=1e3, unit=None,
                        res=256, size=1, format='%3.1f',
                        time_format='%01d ms', proj=False, show=True,
                        show_names=False, title=None):
    """Plot topographic maps of specific time points of evoked data

    Parameters
    ----------
    evoked : Evoked
        The Evoked object.
    times : float | array of floats | None.
        The time point(s) to plot. If None, 10 topographies will be shown
        will a regular time spacing between the first and last time instant.
    ch_type : 'mag' | 'grad' | 'planar1' | 'planar2' | 'eeg'
        The channel type to plot. For 'grad', the gradiometers are collected in
        pairs and the RMS for each pair is plotted.
    layout : None | Layout
        Layout instance specifying sensor positions (does not need to
        be specified for Neuromag data). If possible, the correct layout file
        is inferred from the data; if no appropriate layout file was found, the
        layout is automatically generated from the sensor locations.
    vmin : float | callable
        The value specfying the lower bound of the color range.
        If None, and vmax is None, -vmax is used. Else np.min(data).
        If callable, the output equals vmin(data).
    vmax : float | callable
        The value specfying the upper bound of the color range.
        If None, the maximum absolute value is used. If vmin is None,
        but vmax is not, defaults to np.min(data).
        If callable, the output equals vmax(data).
    cmap : matplotlib colormap
        Colormap. For magnetometers and eeg defaults to 'RdBu_r', else
        'Reds'.
    sensors : bool | str
        Add markers for sensor locations to the plot. Accepts matplotlib plot
        format string (e.g., 'r+' for red plusses).
    colorbar : bool
        Plot a colorbar.
    scale : float | None
        Scale the data for plotting. If None, defaults to 1e6 for eeg, 1e13
        for grad and 1e15 for mag.
    scale_time : float | None
        Scale the time labels. Defaults to 1e3 (ms).
    units : str | None
        The units of the channel types used for colorbar lables. If
        scale == None the unit is automatically determined.
    res : int
        The resolution of the topomap image (n pixels along each side).
    size : float
        Side length per topomap in inches.
    format : str
        String format for colorbar values.
    time_format : str
        String format for topomap values. Defaults to "%01d ms"
    proj : bool | 'interactive'
        If true SSP projections are applied before display. If 'interactive',
        a check box for reversible selection of SSP projection vectors will
        be show.
    show : bool
        Call pyplot.show() at the end.
    show_names : bool | callable
        If True, show channel names on top of the map. If a callable is
        passed, channel names will be formatted using the callable; e.g., to
        delete the prefix 'MEG ' from all channel names, pass the function
        lambda x: x.replace('MEG ', '')
    title : str | None
        Title. If None (default), no title is displayed.
    """
    import matplotlib.pyplot as plt

    if times is None:
        times = np.linspace(evoked.times[0], evoked.times[-1], 10)
    elif np.isscalar(times):
        times = [times]
    if len(times) > 20:
        raise RuntimeError('Too many plots requested. Please pass fewer '
                           'than 20 time instants.')
    tmin, tmax = evoked.times[[0, -1]]
    for t in times:
        if not tmin <= t <= tmax:
            raise ValueError('Times should be between %0.3f and %0.3f. (Got '
                             '%0.3f).' % (tmin, tmax, t))
    if not show_names:
        names = None

    n = len(times)
    nax = n + bool(colorbar)
    width = size * nax
    height = size * 2. + max(0, 0.1 * (4 - size))
    fig = plt.figure(figsize=(width, height))
    w_frame = plt.rcParams['figure.subplot.wspace'] / (2 * nax)
    top_frame = max((0.05 if title is None else 0.15), .2 / size)
    fig.subplots_adjust(left=w_frame, right=1 - w_frame, bottom=0,
                        top=1 - top_frame)
    time_idx = [np.where(evoked.times >= t)[0][0] for t in times]

    vmax_orig, vmin_orig = vmax, vmin
    for jj, key in enumerate(['grad','mag']):

        scale = DEFAULTS['scalings'][key]
        unit = DEFAULTS['units'][key]

        ch_type = key
        picks, pos, merge_grads, names = _prepare_topo_plot(evoked, ch_type,
                                                            layout)
        if proj is True and evoked.proj is not True:
            data = evoked.copy().apply_proj().data
        else:
            data = evoked.data
        data = data[np.ix_(picks, time_idx)] * scale
        
        if merge_grads:
            from mne.layouts.layout import _merge_grad_data
            data = _merge_grad_data(data)

        if type(vmax_orig) is list:
            vmax, vmin = vmax_orig[jj], vmin_orig[jj]
        else:
            vmax, vmin = vmax_orig, vmin_orig
        if vmax is None and vmin is None:
            vmax = np.abs(data).max()
            vmin = -vmax
        else:
            if callable(vmin):
                vmin = vmin(data)
            elif vmin is None:
                vmin = np.min(data)
            if callable(vmax):
                vmax = vmax(data)
            elif vmin is None:
                vmax = np.max(data)

        images = []
        for i, t in enumerate(times):
            plt.subplot(2, nax, (jj-1)*nax + i + 1)
            tp = plot_topomap(data[:, i], pos, vmin=vmin, vmax=vmax, cmap=cmap,
                              sensors=sensors, res=res, names=names,
                              show_names=show_names)
            images.append(tp)
            plt.title(time_format % (t * scale_time))

        if colorbar:
            cax = plt.subplot(2, n + 1, (jj-1)*(n+1) + n + 1)
            plt.colorbar(cax=cax, ticks=[vmin, 0, vmax], format=format)
            # resize the colorbar (by default the color fills the whole axes)
            cpos = cax.get_position()
            cpos.x0 = 1 - (.7 + .1 / size) / nax
            cpos.x1 = cpos.x0 + .1 / nax
            cpos.y0 = jj*0.5 + 0.05
            cpos.y1 = .5*(jj+1) - 0.2
            cax.set_position(cpos)
            if unit is not None:
                cax.set_title(unit)

        if proj == 'interactive':
            _check_delayed_ssp(evoked)
            params = dict(evoked=evoked, fig=fig, projs=evoked.info['projs'],
                          picks=picks, images=images, time_idx=time_idx,
                          scale=scale, merge_grads=merge_grads, res=res, pos=pos,
                          plot_update_proj_callback=_plot_update_evoked_topomap)
            _draw_proj_checkbox(None, params)

    if title is not None:
        plt.suptitle(title, verticalalignment='top', size='x-large')
        tight_layout(pad=2.0, fig=fig)
    if show:
        plt.show()

    return fig

