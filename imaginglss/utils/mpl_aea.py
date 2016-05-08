""" 
    Native matplotlib support of frequently used 2d projections,
    for looking up to the sky.

    This file is initially developed as part of skymapper by Peter Melchior
    based on the example in matplotlib.

    It is later adopted by me (Yu Feng), and I will maintain a copy in
    imaginglss for easier access, also because I do plan to clean up
    the function signatures and variable naming (breaking compatibility with
    old skymapper code).

    The current version adds the ability to generate equal area histograms
    on HealPix pixels.

    It does not depend on healpy, there is a minimal python implementation of 
    healpix at the end of the file; imported in the javascript/lua style.
    
    The intention is one day we will submit a PR of this to matplotlib.

    What does not work:
        
        1. Panning.
        2. Color bar is sometimes in the wrong place
        3. Label locations are poorly calculated.

    What does work:
        Evertying else.

    Author: Yu Feng 
            Peter Melchior

"""
from __future__ import unicode_literals

import matplotlib
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, Polygon
from matplotlib.path import Path
from matplotlib.collections import PolyCollection
from matplotlib.ticker import NullLocator, Formatter, FixedLocator, MaxNLocator
from matplotlib.transforms import Affine2D, BboxTransformTo, Transform, blended_transform_factory, Bbox
from matplotlib.projections import register_projection
import matplotlib.spines as mspines
import matplotlib.axis as maxis

import numpy as np

__author__ = "Yu Feng"
__email__ =  "rainwoodman@gmail.com"

class SkymapperAxes(Axes):
    """
    A base class for a Skymapper axes that takes in ra0, dec0, dec1, dec2.

    The base class takes care of clipping and interpolating with matplotlib.

    Subclass and override class method get_projection_class.

    """
    # The subclass projection must specify a name.  This will be used be the
    # user to select the projection.

    name = None

    @classmethod
    def get_projection_class(kls):
        raise NotImplementedError('Must implement this in subclass')

    def __init__(self, *args, **kwargs):
        self.ra0 = None
        self.dec0 = None
        self.dec1 = None
        self.dec2 = None

        Axes.__init__(self, *args, **kwargs)

        self.cla()

    def _init_axis(self):
        # Axes._init_axis() -- until HammerAxes.xaxis.cla() works.
        self.xaxis = maxis.XAxis(self)
        self.spines['bottom'].register_axis(self.xaxis)
        self.spines['top'].register_axis(self.xaxis)
        self.yaxis = maxis.YAxis(self)
        self.spines['left'].register_axis(self.yaxis)
        self.spines['right'].register_axis(self.yaxis)
        self._update_transScale()

    def cla(self):
        """
        Override to set up some reasonable defaults.
        """
        # Don't forget to call the base class
        Axes.cla(self)

        # Turn off minor ticking altogether
        self.xaxis.set_minor_locator(NullLocator())
        self.yaxis.set_minor_locator(NullLocator())

        self.xaxis.set_major_locator(MaxNLocator(5, prune='both'))
        self.yaxis.set_major_locator(MaxNLocator(5, prune='both'))

        # Do not display ticks -- we only want gridlines and text
        self.xaxis.set_ticks_position('none')
        self.yaxis.set_ticks_position('none')

        self.set_center(None, None)

        # FIXME: probabaly want to override autoscale_view
        # to properly handle wrapping introduced by margin
        # and properlty wrap data. 
        # It doesn't make sense to have xwidth > 360. 
        self._tight = True

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        # There are three important coordinate spaces going on here:
        #
        #    1. Data space: The space of the data itself
        #
        #    2. Axes space: The unit rectangle (0, 0) to (1, 1)
        #       covering the entire plot area.
        #
        #    3. Display space: The coordinates of the resulting image,
        #       often in pixels or dpi/inch.

        # This function makes heavy use of the Transform classes in
        # ``lib/matplotlib/transforms.py.`` For more information, see
        # the inline documentation there.

        # The goal of the first two transformations is to get from the
        # data space (in this case meridian and parallel) to axes
        # space.  It is separated into a non-affine and affine part so
        # that the non-affine part does not have to be recomputed when
        # a simple affine change to the figure has been made (such as
        # resizing the window or changing the dpi).

        # 1) The core transformation from data space into
        # rectilinear space defined in the HammerTransform class.
        self.transProjection = self.get_projection_class()()
        self.transProjection.set_center((180, 0))
        self.transProjection.set_dec1(-65)
        self.transProjection.set_dec2(80)

        # 2) The above has an output range that is not in the unit
        # rectangle, so scale and translate it so it fits correctly
        # within the axes.  The peculiar calculations of xscale and
        # yscale are specific to a Aitoff-Hammer projection, so don't
        # worry about them too much.

        # This will be updated after the xy limits are set.
        self.transAffine = Affine2D()

        # 3) This is the transformation from axes space to display
        # space.
        self.transAxes = BboxTransformTo(self.bbox)

        # Now put these 3 transforms together -- from data all the way
        # to display coordinates.  Using the '+' operator, these
        # transforms will be applied "in order".  The transforms are
        # automatically simplified, if possible, by the underlying
        # transformation framework.
        self.transData = \
            self.transProjection + \
            self.transAffine + \
            self.transAxes

        self.transClip = \
            self.transProjection + \
            self.transAffine

        # The main data transformation is set up.  Now deal with
        # gridlines and tick labels.

        # Longitude gridlines and ticklabels.  The input to these
        # transforms are in display space in x and axes space in y.
        # Therefore, the input values will be in range (-xmin, 0),
        # (xmax, 1).  The goal of these transforms is to go from that
        # space to display space.  The tick labels will be offset 4
        # pixels from the equator.
        self._xaxis_pretransform = \
            Affine2D() \
            .scale(1.0, 180) \
            .translate(0.0, -90)

        self._xaxis_transform = \
            self._xaxis_pretransform + \
            self.transData

        self._xaxis_text1_transform = \
            self._xaxis_pretransform + \
            self.transData + \
            Affine2D().translate(0.0, -8.0)
        self._xaxis_text2_transform = \
            self._xaxis_pretransform+ \
            self.transData + \
            Affine2D().translate(0.0, -8.0)

        # Now set up the transforms for the parallel ticks.  The input to
        # these transforms are in axes space in x and display space in
        # y.  Therefore, the input values will be in range (0, -ymin),
        # (1, ymax).  The goal of these transforms is to go from that
        # space to display space.  The tick labels will be offset 4
        # pixels from the edge of the axes ellipse.
        self._yaxis_stretch = Affine2D().scale(360, 1.0).translate(0.0, 0.0)
        self._yaxis_stretch1 = Affine2D().scale(360, 1.0).translate(0.0, 0.0)
        self._yaxis_stretch2 = Affine2D().scale(360, 1.0).translate(0.0, 0.0)

        self._yaxis_transform = \
            self._yaxis_stretch + \
            self.transData

        self._yaxis_text1_transform = \
            self._yaxis_stretch1 + \
            self.transData
#            Affine2D().translate(-8.0, 0.0)

        self._yaxis_text2_transform = \
            self._yaxis_stretch2 + \
            self.transData
#            Affine2D().translate(8.0, 0.0)

    def _update_affine(self):
        # update the transformations and clip paths
        # after new lims are set.
        if self.ra0 is None:
            x0, x1 = self.viewLim.intervalx
            ra0 = 0.5 * (x0 + x1)
        else:
            ra0 = self.ra0
        if self.dec0 is None:
            y0, y1 = self.viewLim.intervaly
            dec0 = 0.5 * (y0 + y1)
        else:
            dec0 = self.dec0
        if self.dec1 is None:
            y0, y1 = self.viewLim.intervaly
            dec1 = y0 + (y1 - y0) / 12.
        else:
            dec1 = self.dec1
        if self.dec2 is None:
            y0, y1 = self.viewLim.intervaly
            dec2 = y1 - (y1 - y0) / 12.
        else:
            dec2 = self.dec2

        self.transProjection.set_center((ra0, dec0))
        self.transProjection.set_dec1(dec1)
        self.transProjection.set_dec2(dec2)

        self._yaxis_stretch\
            .clear() \
            .scale(self.viewLim.width, 1.0) \
            .translate(self.viewLim.x0, 0)

        self._yaxis_stretch1\
            .clear() \
            .scale(self.viewLim.width, 1.0) \
            .translate(self.viewLim.x0 - 0.00 * self.viewLim.width, 0)

        self._yaxis_stretch2\
            .clear() \
            .scale(self.viewLim.width, 1.0) \
            .translate(self.viewLim.x0 + 0.00 * self.viewLim.width, 0)

        self._xaxis_pretransform \
            .clear() \
            .scale(1.0, self.viewLim.height) \
            .translate(0.0, self.viewLim.y0)

        corners_data = np.array([[self.viewLim.x0, self.viewLim.y0],
                      [ra0,            self.viewLim.y0],
                      [self.viewLim.x1, self.viewLim.y0],
                      [self.viewLim.x1, self.viewLim.y1],
                      [self.viewLim.x0, self.viewLim.y1],])

        corners = self.transProjection.transform_non_affine(corners_data)

        x0 = corners[0][0]
        x1 = corners[2][0]

        # special case when x1 is wrapped back to x0
        # FIXME: I don't think we need it anymore.
        if x0 == x1: x1 = - x0

        y0 = corners[1][1]
        y1 = max([corners[3][1], corners[4][1]])

        xscale = x1 - x0
        yscale = y1 - y0

        self.transAffine.clear() \
            .translate( - (x0 + x1) * 0.5, - (y0 + y1) * 0.5) \
            .scale(0.95 / xscale, 0.95 / yscale)  \
            .translate(0.5, 0.5)

        # now update the clipping path
        path = Path(corners_data)
        path0 = self.transProjection.transform_path(path)
        path = self.transClip.transform_path(path)
        self.patch.set_xy(path.vertices)

    def get_xaxis_transform(self, which='grid'):
        """
        Override this method to provide a transformation for the
        x-axis grid and ticks.
        """
        assert which in ['tick1', 'tick2', 'grid']
        return self._xaxis_transform

    def get_xaxis_text1_transform(self, pixelPad):
        """
        Override this method to provide a transformation for the
        x-axis tick labels.

        Returns a tuple of the form (transform, valign, halign)
        """
        return self._xaxis_text1_transform, 'center', 'center'

    def get_xaxis_text2_transform(self, pixelPad):
        """
        Override this method to provide a transformation for the
        secondary x-axis tick labels.

        Returns a tuple of the form (transform, valign, halign)
        """
        return self._xaxis_text2_transform, 'center', 'center'

    def get_yaxis_transform(self, which='grid'):
        """
        Override this method to provide a transformation for the
        y-axis grid and ticks.
        """
        assert which in ['tick1', 'tick2', 'grid']
        return self._yaxis_transform

    def get_yaxis_text1_transform(self, pixelPad):
        """
        Override this method to provide a transformation for the
        y-axis tick labels.

        Returns a tuple of the form (transform, valign, halign)
        """
        return self._yaxis_text1_transform, 'center', 'center'

    def get_yaxis_text2_transform(self, pixelPad):
        """
        Override this method to provide a transformation for the
        secondary y-axis tick labels.

        Returns a tuple of the form (transform, valign, halign)
        """
        return self._yaxis_text2_transform, 'center', 'center'

    def _gen_axes_patch(self):
        """
        ClipPath.

        Initially set to a size of 2 box in transAxes.

        After xlim and ylim are set, this will be changed to the actual
        region in transData.

        For unclear reason the very initial clip path is always applied
        to the grid. Therefore we set size to 2.0 to avoid bad clipping.
        """
        return Polygon([(0, 0), (2, 0), (2, 2), (0, 2)], fill=False)

    def _gen_axes_spines(self):
        d = {
            'left': mspines.Spine.linear_spine(self, spine_type='left'),
            'right': mspines.Spine.linear_spine(self, spine_type='right'),
            'top': mspines.Spine.linear_spine(self, spine_type='top'),
            'bottom': mspines.Spine.linear_spine(self, spine_type='bottom'),
        }
        d['left'].set_position(('axes', 0))
        d['right'].set_position(('axes', 1))
        d['top'].set_position(('axes', 0))
        d['bottom'].set_position(('axes', 1))
        #FIXME: these spines can be moved wit set_position(('axes', ?)) but
        # 'data' fails. Because the transformation is non-separatable,
        # and because spines / data makes that assumption, we probably
        # do not have a easy way to support moving spines via native matplotlib
        # api on data axis.

        # also the labels currently do not follow the spines. Likely because
        # they are not registered?

        return d

    # Prevent the user from applying scales to one or both of the
    # axes.  In this particular case, scaling the axes wouldn't make
    # sense, so we don't allow it.
    def set_xscale(self, *args, **kwargs):
        if args[0] != 'linear':
            raise NotImplementedError
        Axes.set_xscale(self, *args, **kwargs)

    def set_yscale(self, *args, **kwargs):
        if args[0] != 'linear':
            raise NotImplementedError
        Axes.set_yscale(self, *args, **kwargs)

    def set_center(self, ra0, dec0):
        """ Set the center of ra """
        self.ra0 = ra0
        self.dec0 = dec0
        self._update_affine()

    def set_parallels(self, dec1, dec2):
        """ Set the parallels """
        self.dec1 = dec1
        self.dec2 = dec2
        self._update_affine()

    # when xlim and ylim are updated, the transformation
    # needs to be updated too.
    def set_xlim(self, *args, **kwargs):
        Axes.set_xlim(self, *args, **kwargs)

        # FIXME: wrap x0 x1 to ensure they enclose ra0.
        x0, x1 = self.viewLim.intervalx
        if self.ra0 is not None:
            if not x0 <= self.transProjection.ra0 or \
               not x1 > self.transProjection.ra0:
                raise ValueError("The given limit in RA does not enclose ra0")

        self._update_affine()

    def set_ylim(self, *args, **kwargs):
        Axes.set_ylim(self, *args, **kwargs)
        self._update_affine()

    def histmap(self, ra, dec, weights=None, nside=32, perarea=False, mean=False, **kwargs):
        r = histogrammap(ra, dec, weights, nside, perarea=perarea)

        if weights is not None:
            w, N = r
        else:
            w = r
        if mean:
            mask = N != 0
            w[mask] /= N[mask]
        else:
            mask = w > 0
        return w, mask, self.mapshow(w, mask, nest=False, **kwargs)

    def mapshow(self, map, mask=None, nest=False, **kwargs):
        """ Display a healpix map """
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        defaults = dict(rasterized=True,
                    alpha=0.8,
                    linewidth=0)
        defaults.update(kwargs)
        if mask is None:
            mask = map == map

        if vmin is None:
            vmin = np.nanmin(map[mask])
        if vmax is None:
            vmax = np.nanmax(map[mask])

        coll = HealpixCollection(map, mask, 
                transform=self.transData, **defaults)
        coll.set_clim(vmin=vmin, vmax=vmax)
        self.add_collection(coll)
        self._sci(coll)
        self.autoscale_view(tight=True)
        return coll

    def format_coord(self, lon, lat):
        """
        Override this method to change how the values are displayed in
        the status bar.

        In this case, we want them to be displayed in degrees N/S/E/W.
        """
        lon = lon
        lat = lat
        if lat >= 0.0:
            ns = 'N'
        else:
            ns = 'S'
        if lon >= 0.0:
            ew = 'E'
        else:
            ew = 'W'
        # \u00b0 : degree symbol
        return '%f\u00b0%s, %f\u00b0%s' % (abs(lat), ns, abs(lon), ew)

    class DegreeFormatter(Formatter):
        """
        This is a custom formatter that converts the native unit of
        radians into (truncated) degrees and adds a degree symbol.
        """

        def __init__(self, round_to=1.0):
            self._round_to = round_to

        def __call__(self, x, pos=None):
            degrees = round(x / self._round_to) * self._round_to
            # \u00b0 : degree symbol
            return "%d\u00b0" % degrees

    def set_meridian_grid(self, degrees):
        """
        Set the number of degrees between each meridian grid.

        It provides a more convenient interface to set the ticking than set_xticks would.
        """
        # Set up a FixedLocator at each of the points, evenly spaced
        # by degrees.
        x0, x1 = self.get_xlim()
        number = abs((x1 - x0) / degrees) + 1
        self.xaxis.set_major_locator(
            FixedLocator(
                np.linspace(x0, x1, number, True)[1:-1]))
        # Set the formatter to display the tick labels in degrees,
        # rather than radians.
        self.xaxis.set_major_formatter(self.DegreeFormatter(degrees))

    def set_parallel_grid(self, degrees):
        """
        Set the number of degrees between each meridian grid.

        It provides a more convenient interface than set_yticks would.
        """
        # Set up a FixedLocator at each of the points, evenly spaced
        # by degrees.
        y0, y1 = self.get_ylim()
        number = ((y1 - y0) / degrees) + 1
        self.yaxis.set_major_locator(
            FixedLocator(
                np.linspace(y0, y1, number, True)[1:-1]))
        # Set the formatter to display the tick labels in degrees,
        # rather than radians.
        self.yaxis.set_major_formatter(self.DegreeFormatter(degrees))

    # Interactive panning and zooming is not supported with this projection,
    # so we override all of the following methods to disable it.
    def _in_axes(self, mouseevent):
        if hasattr(self._pan_trans):
            return True
        else:
            return Axes._in_axes(self, mouseevent)

    def can_zoom(self):
        """
        Return True if this axes support the zoom box
        """
        return True

    def start_pan(self, x, y, button):
        self._pan_trans = self.transAxes.inverted() + \
                blended_transform_factory(
                        self._yaxis_stretch,
                        self._xaxis_pretransform,)

    def end_pan(self):
        delattr(self, '_pan_trans')

    def drag_pan(self, button, key, x, y):
        pan1 = self._pan_trans.transform([(x, y)])[0]
        self.set_ra0(360 - pan1[0])
        self.set_dec0(pan1[1])
        self._update_affine()

# now define the Albers equal area axes

class AlbersEqualAreaAxes(SkymapperAxes):
    """
    A custom class for the Albers Equal Area projection.

    https://en.wikipedia.org/wiki/Albers_projection
    """

    name = 'aea'

    @classmethod
    def get_projection_class(kls):
        return kls.AlbersEqualAreaTransform

    # Now, the transforms themselves.
    class AlbersEqualAreaTransform(Transform):
        """
        The base Hammer transform.
        """
        input_dims = 2
        output_dims = 2
        is_separable = False

        def __init__(self, **kwargs):
            Transform.__init__(self, **kwargs)
            self.dec0 = 0
            self.ra0 = 180
            self.dec1 = -60
            self.dec2 = 30
            self._update()

        def set_center(self, center):
            ra0, dec0 = center
            self.ra0  = ra0
            self.dec0 = dec0
            self._update()

        def set_dec1(self, dec1):
            self.dec1 = dec1
            self._update()

        def set_dec2(self, dec2):
            self.dec2 = dec2
            self._update()

        def _update(self):
            self.n = 0.5 * (np.sin(np.radians(self.dec1)) 
                          + np.sin(np.radians(self.dec2)))

            self.C = np.cos(np.radians(self.dec1))**2 + 2 * self.n * np.sin(np.radians(self.dec1))
            self.rho0 = self.__rho__(self.dec0)

        def __rho__(self, dec):
            if self.n == 0:
                return np.sqrt(self.C - 2 * self.n * np.sin(np.radians(dec)))
            else:
                return np.sqrt(self.C - 2 * self.n * np.sin(np.radians(dec))) / self.n

        def transform_non_affine(self, ll):
            """
            Override the transform_non_affine method to implement the custom
            transform.

            The input and output are Nx2 numpy arrays.
            """
            ra = ll[:,0]
            dec = ll[:,1]
            ra0 = self.ra0
            ra_ = np.radians(ra - ra0) # Do not inverse for RA

            # FIXME: problem with the slices sphere: outer parallel needs to be dubplicated at the expense of the central one
            if self.n == 0:
                rt = np.array([
                    self.rho0 * (ra_),
                    - self.rho0 * (np.sin(np.radians(self.dec0) - np.sin(np.radians(dec)))),
                    ]).T
            else:
                theta = self.n * ra_
                rho = self.__rho__(dec)
                rt = np.array([
                       rho*np.sin(theta),
                       self.rho0 - rho*np.cos(theta)]).T
            #if np.isnan(rt).any(): 
            #    raise ValueError('nan occured : ll =%s' % (str(ll)))
            return rt

        # This is where things get interesting.  With this projection,
        # straight lines in data space become curves in display space.
        # This is done by interpolating new values between the input
        # values of the data.  Since ``transform`` must not return a
        # differently-sized array, any transform that requires
        # changing the length of the data array must happen within
        # ``transform_path``.
        def transform_path_non_affine(self, path):
            # Adaptive interpolation:
            # we keep adding control points, till all control points
            # have an error of less than 0.01 (about 1%)
            # or if the number of control points is > 80.
            ra0 = self.ra0
            path = path.cleaned(curves=False)
            v = path.vertices
            diff = v[:, 0] - v[0, 0]
            v00 = v[0][0] - ra0
            while v00 > 180: v00 -= 360
            while v00 < -180: v00 += 360
            v00 += ra0
            v[:, 0] = v00 + diff
            nonstop = path.codes > 0
            path = Path(v[nonstop], path.codes[nonstop])
            isteps = int(path._interpolation_steps * 1.5)
            if isteps < 10: isteps = 10
            while True:
                ipath = path.interpolated(isteps)
                tiv = self.transform(ipath.vertices)
                itv = Path(self.transform(path.vertices)).interpolated(isteps).vertices
                if np.mean(np.abs(tiv - itv)) < 0.01:
                    break
                if isteps > 20:
                    break
                isteps = int(isteps * 1.5)
            return Path(tiv, ipath.codes)

        transform_path_non_affine.__doc__ = \
            Transform.transform_path_non_affine.__doc__

        if matplotlib.__version__ < '1.2':
            transform = transform_non_affine
            transform_path = transform_path_non_affine
            transform_path.__doc__ = Transform.transform_path.__doc__

        def inverted(self):
            return AlbersEqualAreaAxes.InvertedAlbersEqualAreaTransform(self)
        inverted.__doc__ = Transform.inverted.__doc__

    class InvertedAlbersEqualAreaTransform(Transform):
        """ Inverted transform.

            This will always only give values in the prime ra0-180 ~ ra0+180 range, I believe.
            So it is inherently broken. I wonder when matplotlib actually calls this function,
            given that interactive is disabled.
        """
        input_dims = 2
        output_dims = 2
        is_separable = False

        def __init__(self, inverted, **kwargs):
            Transform.__init__(self, **kwargs)
            self.inverted = inverted

        def transform_non_affine(self, xy):
            x = xy[:,0]
            y = xy[:,1]
            inverted = self.inverted

            rho = np.sqrt(x**2 + (inverted.rho0 - y)**2)

            # make sure that the signs are correct
            if inverted.n == 0:
                rt = np.degrees(
                        [
                    np.radians(inverted.ra0) + x / inverted.rho0,
                    np.arcsin(y / inverted.rho0 + np.sin(np.radians(inverted.dec0)))
                        ]).T
                return rt
            elif inverted.n > 0:
                theta = np.degrees(np.arctan2(x, inverted.rho0 - y))
            else:
                theta = np.degrees(np.arctan2(-x, -(inverted.rho0 - y)))
            return np.degrees([np.radians(inverted.ra0) + theta/inverted.n,
                np.arcsin((inverted.C - (rho * inverted.n)**2)/(2*inverted.n))]).T

            transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

        if matplotlib.__version__ < '1.2':
            transform = transform_non_affine

        def inverted(self):
            # The inverse of the inverse is the original transform... ;)
            return self.inverted

        inverted.__doc__ = Transform.inverted.__doc__

class HealpixCollection(PolyCollection):
    def __init__(self, map, mask, nest=False, **kwargs):
        self.v = _boundary(mask, nest)
        PolyCollection.__init__(self, self.v, array=map[mask], **kwargs)

    def get_datalim(self, transData):
        """ The data lim of a healpix collection.
        """ 
        # FIXME: it is currently set to the full sky.
        #    This could have been trimmed down. 
        #    We want to set xlim smartly such that the largest
        #    empty region is chopped off. I think it is possible, by
        #    doing a histogram in ra, for example. 
        vmin = (0, -90)
        vmax = (360, 90)
        return Bbox((vmin, vmax))

# a few helper functions talking to healpy/healpix.
def _boundary(mask, nest=False):
    """Generate healpix vertices for pixels where mask is True

    Args:
        pix: list of pixel numbers
        nest: nested or not
        nside: HealPix nside

    Returns:
        vertices
        vertices: (N,4,2), RA/Dec coordinates of 4 boundary points of cell
    """

    pix = mask.nonzero()[0]

    nside = healpix.npix2nside(len(mask))

    vertices = np.zeros((pix.size, 4, 2))
    theta, phi = healpix.vertices(nside, pix)
    theta = np.degrees(theta)
    phi = np.degrees(phi)
    diff = phi - phi[:, 0][:, None]
    diff[diff > 180] -= 360
    diff[diff < -180] += 360
    phi = phi[:, 0][:, None] + diff
    vertices[:,:,0] = phi
    vertices[:,:,1] = 90.0 - theta

    return vertices

def histogrammap(ra, dec, weights=None, nside=32, perarea=False):
    ipix = healpix.ang2pix(nside, np.radians(90-dec), np.radians(ra))
    npix = healpix.nside2npix(nside)
    if perarea:
        npix = healpix.nside2npix(nside)
        sky = 360. ** 2 / np.pi
        area = 1. * (sky / npix)
    else:
        area = 1

    if weights is not None:
        w = np.bincount(ipix, weights=weights, minlength=npix)
        N = np.bincount(ipix, minlength=npix)
        w = w / area
        N = N / area
        return w, N
    else:
        w = 1.0 * np.bincount(ipix, minlength=npix)
        return w / area

# Now register the projection with matplotlib so the user can select
# it.
register_projection(AlbersEqualAreaAxes)

def create_healpix():
    """ A pure python (numpy-based) version of key healpix functions.

        The ring scheme is implemented. 

        Depencency: numpy.

        It shall probably be self-hosted as an individual python package.

        Author: Yu Feng <rainwoodman@gmail.com>
    """

    import numpy

    def npix2nside(npix):
        # FIXME: this could be buggy for large npix
        nside2 = npix // 12
        nside = numpy.array(nside2 ** 0.5).astype('i8')
        return nside

    def nside2npix(nside):
        return nside * nside * 12

    def ang2pix(nside, theta, phi):
        r"""Convert angle :math:`\theta` :math:`\phi` to pixel.

            This is translated from chealpix.c; but refer to Section 4.1 of
            http://adsabs.harvard.edu/abs/2005ApJ...622..759G
        """
        nside, theta, phi = numpy.lib.stride_tricks.broadcast_arrays(nside, theta, phi)
        
        def equatorial(nside, tt, z):
            t1 = nside * (0.5 + tt)
            t2 = nside * z * 0.75
            jp = (t1 - t2).astype('i8')
            jm = (t1 + t2).astype('i8')
            ir = nside + 1 + jp - jm # in {1, 2n + 1}
            kshift = 1 - (ir & 1) # kshift=1 if ir even, 0 odd 
     
            ip = (jp + jm - nside + kshift + 1) // 2 # in {0, 4n - 1}
            
            ip = ip % (4 * nside)
            return nside * (nside - 1) * 2 + (ir - 1) * 4 * nside + ip
            
        def polecaps(nside, tt, z, s):
            tp = tt - numpy.floor(tt)
            za = numpy.abs(z)
            tmp = nside * s / ((1 + za) / 3) ** 0.5
            mp = za > 0.99
            tmp[mp] = nside[mp] * (3 *(1-za[mp])) ** 0.5
            jp = (tp * tmp).astype('i8')
            jm = ((1 - tp) * tmp).astype('i8')
            ir = jp + jm + 1
            ip = (tt * ir).astype('i8')
            ip = ip % (4 * ir)

            r1 = 2 * ir * (ir - 1) 
            r2 = 2 * ir * (ir + 1)
     
            r = numpy.empty_like(r1)
            
            r[z > 0] = r1[z > 0] + ip[z > 0]
            r[z < 0] = 12 * nside[z < 0] * nside[z < 0] - r2[z < 0] + ip[z < 0]
            return r
        
        z = numpy.cos(theta)
        s = numpy.sin(theta)
        
        tt = (phi / (0.5 * numpy.pi) ) % 4 # in [0, 4]
        
        result = numpy.zeros(z.shape, dtype='i8')
        mask = (z < 2. / 3) & (z > -2. / 3)
      
        result[mask] = equatorial(nside[mask], tt[mask], z[mask])
        result[~mask] = polecaps(nside[~mask], tt[~mask], z[~mask], s[~mask])
        return result
        
    def pix2ang(nside, pix):
        r"""Convert pixel to angle :math:`\theta` :math:`\phi`.

            nside and pix are broadcast with numpy rules.

            Returns: theta, phi

            This is translated from chealpix.c; but refer to Section 4.1 of
            http://adsabs.harvard.edu/abs/2005ApJ...622..759G
        """
        nside, pix = numpy.lib.stride_tricks.broadcast_arrays(nside, pix)
        
        ncap = nside * (nside - 1) * 2
        npix = 12 * nside * nside
        
        def northpole(pix, npix):
            iring = (1 + ((1 + 2 * pix) ** 0.5)).astype('i8') // 2
            iphi = (pix + 1) - 2 * iring * (iring - 1)
            z = 1.0 - (iring*iring) * 4. / npix
            phi = (iphi - 0.5) * 0.5 * numpy.pi / iring
            return z, phi
        
        def equatorial(pix, nside, npix, ncap):
            ip = pix - ncap
            iring = ip // (4 * nside) + nside
            iphi = ip % (4 * nside) + 1
            fodd = (((iring + nside) &1) + 1.) * 0.5
            z = (2 * nside - iring) * nside * 8.0 / npix
            phi = (iphi - fodd) * (0.5 * numpy.pi) / nside
            return z, phi
        
        def southpole(pix, npix):
            ip = npix - pix
            iring = (1 + ((2 * ip - 1)**0.5).astype('i8')) // 2
            iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1))
            z = -1 + (iring * iring) * 4. / npix
            phi = (iphi - 0.5 ) * 0.5 * numpy.pi / iring
            return z, phi
        
        mask1 = pix < ncap
        
        mask2 = (~mask1) & (pix < npix - ncap)
        mask3 = pix >= npix - ncap

        z = numpy.zeros(pix.shape, dtype='f8')
        phi = numpy.zeros(pix.shape, dtype='f8')
        z[mask1], phi[mask1] = northpole(pix[mask1], npix[mask1])
        z[mask2], phi[mask2] = equatorial(pix[mask2], nside[mask2], npix[mask2], ncap[mask2])
        z[mask3], phi[mask3] = southpole(pix[mask3], npix[mask3])
        return numpy.arccos(z), phi

    def ang2xy(theta, phi):
        r"""Convert :math:`\theta` :math:`\phi` to :math:`x_s` :math:`y_s`.

            Returns: x, y

            Refer to Section 4.4 of http://adsabs.harvard.edu/abs/2005ApJ...622..759G
        """
        theta, phi = numpy.lib.stride_tricks.broadcast_arrays(theta, phi)
        z = numpy.cos(theta)
        x = numpy.empty(theta.shape, dtype='f8')
        y = numpy.empty(theta.shape, dtype='f8')
        def sigma(z):
            return numpy.sign(z) * (2 - (3 * (1- numpy.abs(z))) ** 0.5)
                
        def equatorial(z, phi):
            return phi, 3 * numpy.pi / 8 * z
        def polarcaps(z, phi):
            s = sigma(z)
            x = phi - (numpy.abs(s) - 1) * (phi % (0.5 * numpy.pi) - 0.25 * numpy.pi)
            y = 0.25 * numpy.pi * s
            return x, y
        
        mask = numpy.abs(z) < 2. / 3

        x[mask], y[mask] = equatorial(z[mask], phi[mask])
        x[~mask], y[~mask] = polarcaps(z[~mask], phi[~mask])
        return x, y

    def xy2ang(x, y):
        r"""Convert :math:`x_s` :math:`y_s` to :math:`\theta` :math:`\phi`.
            
            Returns: theta, phi

            Refer to Section 4.4 of http://adsabs.harvard.edu/abs/2005ApJ...622..759G
        """
        x, y = numpy.lib.stride_tricks.broadcast_arrays(x, y)
        
        theta = numpy.empty(x.shape, dtype='f8')
        phi = numpy.empty(x.shape, dtype='f8')
        
        def equatorial(x, y):
            return numpy.arccos(8 * y / (3 * numpy.pi)), x
        
        def polarcaps(x, y):
            ya = numpy.abs(y)
            xt = x % (0.5 * numpy.pi)
            phi = x - (ya - numpy.pi * 0.25) / (ya - numpy.pi * 0.5) * (xt - 0.25 * numpy.pi)
            z = (1 - 1. / 3 * (2 - 4 * ya / numpy.pi)**2) * y / ya
            return numpy.arccos(z), phi
        
        mask = numpy.abs(y) < numpy.pi * 0.25
       
        theta[mask], phi[mask] = equatorial(x[mask], y[mask])
        theta[~mask], phi[~mask] = polarcaps(x[~mask], y[~mask])
        return theta, phi

    def vertices(nside, pix):
        r""" Calculate the vertices for pixels 

            Returns: theta, phi
                for each (nside, pix) pair, a four-vector of theta, and
                a four-vector of phi is returned, corresponding to
                the theta, phi of each vertex of the pixel boundary.
        """
        nside, pix = numpy.lib.stride_tricks.broadcast_arrays(nside, pix)
        x = numpy.zeros(nside.shape, dtype=('f8', 4))
        y = numpy.zeros(nside.shape, dtype=('f8', 4))
        theta, phi = pix2ang(nside, pix)
        xc, yc = ang2xy(theta, phi)
        xstep = numpy.pi / (2 * nside)
        ystep = numpy.pi / (2 * nside)
        x[..., 0] = xc - 0.5 * xstep
        y[..., 0] = yc
        x[..., 1] = xc
        y[..., 1] = yc + 0.5 * ystep
        x[..., 2] = xc + 0.5 * xstep
        y[..., 2] = yc
        x[..., 3] = xc
        y[..., 3] = yc - 0.5 * ystep
        
        theta, phi = xy2ang(x, y)
        return theta, phi
    return locals()

class Namespace(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

healpix = Namespace(**create_healpix())

if __name__ == '__main__':
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    # Now make a simple example using the custom projection.

    import numpy as np

    fig = Figure(figsize=(6, 6))

#    ra = np.random.uniform(size=100, low=0, high=360)
#    dec = np.random.uniform(size=100, low=-90, high=90)
    ra = np.linspace(0, 360, 100)
    dec = np.linspace(-90, 90, 100)
    ax = fig.add_subplot(111, projection="aea")
    ax.set_xlim(359, 0)
    ax.set_ylim(-70, 70)
    ax.set_parallels(-20, 60)
    ax.set_center(180, 0)
    ax.plot(ra, dec, '*')
    ax.axhline(-20)
    ax.axvline(140)

    ra = np.random.uniform(size=1000, low=0, high=360)
    dec = np.random.uniform(size=1000, low=-90, high=90)
    ax.histmap(ra, dec, nside=8)

    ax.tick_params(labelright=True, labeltop=True)

    ax.set_meridian_grid(30)
    ax.set_parallel_grid(30)
    ax.grid()
    ax.tricontour(ra, dec, ra)
    fig.colorbar(ax._gci())
    canvas = FigureCanvasAgg(fig)
    fig.savefig('xxx.png')
