# ChromaPlot
Chromaticity Diagram Plotter for LaTex

Previous attempts:

[Paolo Brasolin](https://paolobrasolin.github.io/) first published his attempt at a `pgfplots` based vectorized chromaticity diagram [on Stackexchange](https://tex.stackexchange.com/questions/177079/tikz-chromaticity-diagram).

### Python Packages

| Package  | Author | Relevant Namespace/Module | Documentation | Source |
| -------- | -------- | -------- |-------- | ------- |
| LuxPy | Kevin Smet | [`luxpy.color.utils.plotters`](https://github.com/ksmet1977/luxpy/blob/master/luxpy/color/utils/plotters.py) | [readthedocs](https://ksmet1977.github.io/luxpy/build/html/index.html) | [Github](https://github.com/ksmet1977/luxpy) |
| Colour | Various Authors| [`colour.plotting`](https://github.com/colour-science/colour/tree/develop/colour/plotting) |  [readthedocs](https://colour.readthedocs.io/en/latest/) | [Github](https://github.com/colour-science/colour) |
| ColorPy | Mark Kness | ??? |  [Website](http://markkness.net/colorpy/ColorPy.html) | [Github](https://github.com/markkness/ColorPy) |
| colorio | Nico Schlömer | ??? |  [PyPi](https://pypi.org/project/colorio/) | [Github](https://github.com/nschloe/colorio/) |
| python-colormath | Gregory Taylor | ??? |  [readthedocs](https://python-colormath.readthedocs.io/en/latest/) | [Github](https://github.com/gtaylor/python-colormath) |


#### Chromaticity diagram plotting in LuxPy

Plotting functions are defined in [`luxpy/blob/master/luxpy/color/utils/plotters.py`](https://github.com/ksmet1977/luxpy/blob/master/luxpy/color/utils/plotters.py).

The chromaticity diagram is plotted by the function []`luxpy.color.utils.plot_chromaticity_diagram_colors`](https://ksmet1977.github.io/luxpy/build/html/color.html?highlight=plot_chromaticity_diagram_colors#luxpy.color.utils.plot_chromaticity_diagram_colors).

```
def plot_chromaticity_diagram_colors(diagram_samples = 256, diagram_opacity = 1.0, diagram_lightness = 0.25,\
                                      cieobs = _CIEOBS, cspace = 'Yxy', cspace_pars = {},\
                                      show = True, axh = None,\
                                      show_grid = True, label_fontname = 'Times New Roman', label_fontsize = 12,\
                                      **kwargs):
    """
    Plot the chromaticity diagram colors.

    Args:
        :diagram_samples:
            | 256, optional
            | Sampling resolution of color space.
        :diagram_opacity:
            | 1.0, optional
            | Sets opacity of chromaticity diagram
        :diagram_lightness:
            | 0.25, optional
            | Sets lightness of chromaticity diagram
        :axh:
            | None or axes handle, optional
            | Determines axes to plot data in.
            | None: make new figure.
        :show:
            | True or False, optional
            | Invoke matplotlib.pyplot.show() right after plotting
        :cieobs:
            | luxpy._CIEOBS or str, optional
            | Determines CMF set to calculate spectrum locus or other.
        :cspace:
            | luxpy._CSPACE or str, optional
            | Determines color space / chromaticity diagram to plot data in.
            | Note that data is expected to be in specified :cspace:
        :cspace_pars:
            | {} or dict, optional
            | Dict with parameters required by color space specified in :cspace:
            | (for use with luxpy.colortf())
        :show_grid:
            | True, optional
            | Show grid (True) or not (False)
        :label_fontname:
            | 'Times New Roman', optional
            | Sets font type of axis labels.
        :label_fontsize:
            | 12, optional
            | Sets font size of axis labels.
        :kwargs:
            | additional keyword arguments for use with matplotlib.pyplot.

    Returns:

    """
    offset = _EPS
    ii, jj = np.meshgrid(np.linspace(offset, 1 + offset, diagram_samples), np.linspace(1+offset, offset, diagram_samples))
    ij = np.dstack((ii, jj))

    SL =  _CMF[cieobs]['bar'][1:4].T
    SL = np.vstack((SL,SL[0]))
    SL = 100.0*SL/SL[:,1,None]
    SL = colortf(SL, tf = cspace, tfa0 = cspace_pars)
    Y,x,y = asplit(SL)
    SL = np.vstack((x,y)).T


    ij2D = ij.reshape((diagram_samples**2,2))
    ij2D = np.hstack((diagram_lightness*100*np.ones((ij2D.shape[0],1)), ij2D))
    xyz = colortf(ij2D, tf = cspace + '>xyz', tfa0 = cspace_pars)

    xyz[xyz < 0] = 0
    xyz[np.isinf(xyz.sum(axis=1)),:] = np.nan
    xyz[np.isnan(xyz.sum(axis=1)),:] = offset

    srgb = xyz_to_srgb(xyz)
    srgb = srgb/srgb.max()
    srgb = srgb.reshape((diagram_samples,diagram_samples,3))

    if show == True:
        if axh is None:
            fig = plt.figure()
            axh = fig.add_subplot(111)
        polygon = Polygon(SL, facecolor='none', edgecolor='none')
        axh.add_patch(polygon)
        image = axh.imshow(
            srgb,
            interpolation='bilinear',
            extent = (0.0, 1, -0.05, 1),
            clip_path=None,
            alpha=diagram_opacity)
        image.set_clip_path(polygon)
        plt.plot(x,y, color = 'darkgray')
        if cspace == 'Yxy':
            plt.xlim([0,1])
            plt.ylim([0,1])
        elif cspace == 'Yuv':
            plt.xlim([0,0.6])
            plt.ylim([0,0.6])
        if (cspace is not None):
            xlabel = _CSPACE_AXES[cspace][1]
            ylabel = _CSPACE_AXES[cspace][2]
            if (label_fontname is not None) & (label_fontsize is not None):
                plt.xlabel(xlabel, fontname = label_fontname, fontsize = label_fontsize)
                plt.ylabel(ylabel, fontname = label_fontname, fontsize = label_fontsize)

        if show_grid == True:
            plt.grid()
        #plt.show()

        return axh
    else:
        return None
```


#### Chromaticity diagram plotting in Colour

Plotting functions are defined in [`colour/plotting/diagrams.py`](https://github.com/colour-science/colour/blob/develop/colour/plotting/diagrams.py).

The chromaticity diagram is loaded as background into a `matplotlib` figure from file: `../_static/Plotting_Plot_Chromaticity_Diagram_CIE1931.png`







### Theory

#### Chromaticity Diagrams

- [Explanation of the "star" in the chromaticity diagram in the description of a Wikimedia file](https://commons.wikimedia.org/wiki/File:CIExy1931_srgb_gamut.png)
-

### Glossary

[IES Definitions](https://www.ies.org/definitions/spectrum-locus/)

| Terms  | Definition | Source | Additional Source |
| -------- | -------- | ------------- | ---------- |
| Spectral [locus](https://en.wikipedia.org/wiki/Locus_(mathematics)) | The locus of points representing the colors of the visible spectrum in a chromaticity diagram. | [IES](https://www.ies.org/definitions/spectrum-locus/) | |
| Saturation | Test | Some CIE Source | [Wikipedia Article](https://en.wikipedia.org/wiki/Colorfulness) |
| Luminance | Test | Some CIE Source | [Wikipedia Article](https://en.wikipedia.org/wiki/Luminance)|



### Datasets

[UCL Color Vision and Research Laboratory](http://www.cvrl.org/)
