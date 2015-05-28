import numpy as np
import pandas as pd
import sys
import ast


class GridProc(object):
    """
    Container class for holding and manipulating gridded VIC routing data.

    Attributes
    ----------
    bbox : The geographical bounding box for viewing the gridded data
           (xmin, ymin, xmax, ymax)
    shape : The shape of the gridded data (nrows, ncolumns)
    cellsize : The length/width of each grid cell (assumed to be square).
    grid_props : dict containing metadata for each gridded dataset.

    Methods
    -------
    read_input : add a gridded dataset (dem, flowdir, accumulation)
                 to flow_grid instance.
    nearest_cell : Returns the index (column, row) of the cell closest
                   to a given geographical coordinate (x, y).
    flowdir : Generate a flow direction grid from a given digital elevation
              dataset (dem). Does not currently handle flats.
    catchment : Delineate the watershed for a given pour point (x, y)
                or (column, row).
    fraction : Generate the fractional contributing area for a coarse
               scale flow direction grid based on a fine-scale flow
               direction grid.

    Reserved Datasets
    -----------------
    dem : digital elevation grid
    dir : flow direction grid
    acc : flow accumulation grid
    catch : Catchment delineated from 'dir' and a given pour point
    frac : fractional contributing area grid
    """

    def __init__(self):
        self.grid_props = {}

    def add_data(self, data, data_name, bbox=None, shape=None, cellsize=None,
            crs=None, nodata=None, **kwargs):
        """
        A generic method for adding data into a FlowGrid instance.
        Inserts data into a named attribute of FlowGrid (name of attribute
        determined by 'data_name').

        Parameters
        ----------
        data : numpy ndarray
               If data is from a file, 'input_type' should be set to the
               appropriate value ('ascii' or 'raster').
        data_name : 'dem', 'dir', 'acc' or other string
                     Name of dataset. Will determine the name of the attribute
                     representing the gridded data. Default values are used
                     internally by some class methods:
                         'dem' : digital elevation data
                         'dir' : flow direction data
                         'acc' : flow accumulation (upstream area) data

        bbox : tuple
               Bounding box of data.
        shape : tuple
                Shape (rows, columns) of data.
        cellsize : float or int
                   Cellsize of gridded data.
        crs : dict
              Coordinate reference system of gridded data.
        nodata : int or float
                 Value indicating no data.

        """
        if not isinstance(data, np.ndarray):
            raise TypeError('Input data must be ndarray')

        # if there are no datasets, initialize bbox, shape,
        # cellsize and crs based on incoming data
        if len(self.grid_props) < 1:

            # check validity of bbox
            if ((hasattr(bbox, "__len__")) and (not isinstance(bbox, str))
                    and (len(bbox) == 4)):
                bbox = tuple(bbox)
            else:
                raise TypeError('bbox must be a tuple of length 4.')
            # check validity of shape
            if ((hasattr(shape, "__len__")) and (not isinstance(shape, str))
                    and (len(shape) == 2) and (isinstance(sum(shape), int))):
                shape = tuple(shape)
            else:
                raise TypeError('shape must be a tuple of ints of length 2.')
            # check validity of cellsize
            if not isinstance(cellsize, (int, float)):
                raise TypeError('cellsize must be an int or float.')

            # initialize instance metadata
            self.bbox = bbox
            self.shape = shape
            self.cellsize = cellsize
            self.crs = crs
            self.mask = np.ones(self.shape, dtype=np.bool)
            self.shape_min = np.min_scalar_type(max(self.shape))
            self.size_min = np.min_scalar_type(data.size)

        # if there are existing datasets, conform incoming
        # data to bbox
        else:
            try:
                np.testing.assert_almost_equal(cellsize, self.cellsize)
            except:
                raise AssertionError('Grid cellsize not equal')

        # assign new data to attribute; record nodata value
        self.grid_props.update({data_name : {}})
        self.grid_props[data_name].update({'bbox' : bbox})
        self.grid_props[data_name].update({'shape' : shape})
        self.grid_props[data_name].update({'cellsize' : cellsize})
        self.grid_props[data_name].update({'nodata' : nodata})
        self.grid_props[data_name].update({'crs' : crs})
        setattr(self, data_name, data)

    def read_ascii(self, data, data_name, skiprows=6, **kwargs):
        """
        Reads data from an ascii file into a named attribute of flow_grid
        (name of attribute determined by 'data_name').

        Parameters
        ----------
        data : File name (string) or numpy ndarray
               If data is from a file, 'input_type' should be set to the
               appropriate value ('ascii' or 'raster').
        data_name : 'dem', 'dir', 'acc' or other string
                     Name of dataset. Will determine the name of the attribute
                     representing the gridded data. Default values are used
                     internally by some class methods:
                         'dem' : digital elevation data
                         'dir' : flow direction data
                         'acc' : flow accumulation (upstream area) data
        skiprows : The number of rows taken up by the header.

        Additional keyword arguments are passed to numpy.loadtxt()
        """

        with open(data) as header:
            ncols = int(header.readline().split()[1])
            nrows = int(header.readline().split()[1])
            xll = ast.literal_eval(header.readline().split()[1])
            yll = ast.literal_eval(header.readline().split()[1])
            cellsize = ast.literal_eval(header.readline().split()[1])
            nodata = ast.literal_eval(header.readline().split()[1])
            shape = (nrows, ncols)
            bbox = (xll, yll, xll + ncols * cellsize, yll + nrows * cellsize)
        data = np.loadtxt(data, skiprows=skiprows, **kwargs)
        nodata = data.dtype.type(nodata)
        self.add_data(data, data_name, bbox, shape, cellsize, crs, nodata)

    def read_raster(self, data_name, band=1, **kwargs):
        """
        Reads data from a raster file into a named attribute of flow_grid
        (name of attribute determined by 'data_name').

        Parameters
        ----------
        data : File name (string) or numpy ndarray
               If data is from a file, 'input_type' should be set to the
               appropriate value ('ascii' or 'raster').
        data_name : 'dem', 'dir', 'acc' or other string
                     Name of dataset. Will determine the name of the attribute
                     representing the gridded data. Default values are used
                     internally by some class methods:
                         'dem' : digital elevation data
                         'dir' : flow direction data
                         'acc' : flow accumulation (upstream area) data
        band : int
               The band number to read.

        Additional keyword arguments are passed to rasterio.open()
        """
        # read raster file
        import rasterio
        f = rasterio.open(data, **kwargs)
        crs = f.crs
        bbox = tuple(f.bounds)
        shape = f.shape
        cellsize = f.affine[0]
        nodata = f.nodatavals[0]
        if len(f.indexes) > 1:
            data = np.ma.filled(f.read_band(band))
        else:
            data = np.ma.filled(f.read())
            f.close()
            data = data.reshape(shape)
        nodata = data.dtype.type(nodata)
        self.add_data(data, data_name, bbox, shape, cellsize, crs, nodata)

    def bbox_indices(self, bbox, shape, precision=7):
        """
        Return row and column coordinates of a bounding box at a
        given cellsize.

        Parameters
        ----------
        bbox : tuple
               bbox of new data
        shape : tuple
                the shape of the 2D array (rows, columns)
        precision : int
                    Precision to use when matching geographic coordinates.
        """
        rows = np.around(np.linspace(bbox[1], bbox[3],
               shape[0], endpoint=False)[::-1], precision)
        cols = np.around(np.linspace(bbox[0], bbox[2],
               shape[1], endpoint=False), precision)
        return rows, cols

    def view(self, data_name, mask=True):
        """
        Return a copy of a gridded dataset clipped to the bounding box
        (self.bbox) with cells outside the catchment mask (self.mask)
        displayed as 'nodata' (self.grid_props[data_name]['nodata']).

        Parameters
        ----------
        data_name : string
                    Name of the dataset to be viewed.
        mask : bool
               Whether or not to "mask" the view using self.mask.
        """
        selfrows, selfcols = self.bbox_indices(self.bbox, self.shape)
        rows, cols = self.bbox_indices(self.grid_props[data_name]['bbox'],
                                       self.grid_props[data_name]['shape'])
        outview = (pd.DataFrame(getattr(self, data_name),
                                index=rows, columns=cols)
                   .reindex(selfrows).reindex_axis(selfcols, axis=1)
                   .fillna(self.grid_props[data_name]['nodata']).values)
        if mask:
            return np.where(self.mask,
                            outview,
                            self.grid_props[data_name]['nodata'])
        else:
            return outview

    def nearest_cell(self, x, y):
        """
        Returns the index of the cell (column, row) closest
        to a given geographical coordinate.

        Parameters
        ----------
        x : int or float
            x coordinate.
        y : int or float
            y coordinate.
        """

        # create coordinate grid of cell centroids based on bbox
        coords = np.meshgrid(
            np.linspace(self.bbox[0] + self.cellsize / 2.0,
                        self.bbox[2] + self.cellsize / 2.0,
                        self.shape[1], endpoint=False),
            np.linspace(self.bbox[1] + self.cellsize / 2.0,
                        self.bbox[3] + self.cellsize / 2.0,
                        self.shape[0], endpoint=False)[::-1])

        # select nearest cell based on euclidian distance
        nearest = np.unravel_index(
            np.argmin(np.sqrt((coords[0] - x) ** 2 + (coords[1] - y) ** 2)),
            self.shape)
        return nearest[1], nearest[0]

    def flowdir(self, data_name=None, include_edges=True, nodata=0, flat=-1,
            dirmap=(1, 2, 3, 4, 5, 6, 7, 8), inplace=True):
        """
        Generates a flow direction grid from a DEM grid.

        Parameters
        ----------
        data : numpy ndarray
               Array representing DEM grid
        include_edges : bool
                        Whether to include outer rim of grid.
        nodata : int
                 value to indicate nodata in output array.
        flat : int
               value to indicate flat areas in output array.
        dirmap : list or tuple (length 8)
                 List of integer values representing the following
                 cardinal and intercardinal directions (in order):
                 [N, NE, E, SE, S, SW, W, NW]
        inplace : bool
                  If True, write output array to self.dir
        """

        if len(dirmap) != 8:
            raise AssertionError('dirmap must be a sequence of length 8')

        # if data not provided, use self.dem
        if data_name is None:
            if hasattr(self, 'dem'):
                data = self.view('dem', mask=False)
        elif isinstance(data_name, str):
            data = self.view(data_name, mask=False)
        else:
            raise TypeError('data_name must be string type.')

        # generate grid of indices
        indices = np.indices(self.shape, dtype=self.shape_min)

        # handle nodata values in dem
        dem_nodata = self.grid_props['dem']['nodata']
        dem_mask = (data == dem_nodata)
        np.place(data, dem_mask, np.iinfo(data.dtype.type).max)

        # initialize indices of corners
        corners = {
        'nw' : {'k' : tuple(indices[:, 0, 0]),
                'v' : [[0, 1, 1],  [1, 1, 0]],
                'pad': np.array([3, 4, 5])},
        'ne' : {'k' : tuple(indices[:, 0, -1]),
                'v' : [[1, 1, 0],  [-1, -2, -2]],
                'pad': np.array([5, 6, 7])},
        'sw' : {'k' : tuple(indices[:, -1, 0]),
                'v' : [[-2, -2, -1],  [0, 1, 1]],
                'pad': np.array([1, 2, 3])},
        'se' : {'k' : tuple(indices[:, -1, -1]),
                'v' : [[-1, -2, -2],  [-2, -2, -1]],
                'pad': np.array([7, 8, 1])}
        }

        # initialize indices of edges
        edges = {
        'n' : {'k' : tuple(indices[:, 0, 1:-1]),
               'pad' : np.array([3, 4, 5, 6, 7])},
        'w' : {'k' : tuple(indices[:, 1:-1, 0]),
               'pad' : np.array([1, 2, 3, 4, 5])},
        'e' : {'k' : tuple(indices[:, 1:-1, -1]),
               'pad' : np.array([1, 5, 6, 7, 8])},
        's' : {'k' : tuple(indices[:, -1, 1:-1]),
               'pad' : np.array([1, 2, 3, 7, 8])}
        }

        # initialize indices of body (all cells except edges and corners)
        body = indices[:, 1:-1, 1:-1]

        # initialize output array
        outmap = np.full(self.shape, nodata, dtype=np.int8)

        # for each entry in "body" determine flow direction based
        # on steepest neighboring slope
        for i, j in np.nditer(tuple(body), flags=['external_loop']):
            dat = data[i, j]
            sur = data[self._select_surround(i, j)]
            a = ((dat - sur) > 0).any(axis=0)
            b = np.argmax((dat - sur), axis=0) + 1
            c = flat
            outmap[i, j] = np.where(a, b, c)

        # determine flow direction for edges and corners, if desired
        if include_edges:

            # fill corners
            for corner in corners.keys():
                dat = data[corners[corner]['k']]
                sur = data[corners[corner]['v']]
                if ((dat - sur) > 0).any():
                    outmap[corners[corner]['k']] = \
                            corners[corner]['pad'][np.argmax(dat - sur)]
                else:
                    outmap[corners[corner]['k']] = flat

            # fill edges
            for edge in edges.keys():
                dat = data[edges[edge]['k']]
                sur = data[self._select_edge_sur(edges, edge)]
                a = ((dat - sur) > 0).any(axis=0)
                b = edges[edge]['pad'][np.argmax((dat - sur), axis=0)]
                c = flat
                outmap[edges[edge]['k']] = np.where(a, b, c)

        # If direction numbering isn't default, convert values of output array.
        if dirmap != (1, 2, 3, 4, 5, 6, 7, 8):
            dir_d = dict(zip((1, 2, 3, 4, 5, 6, 7, 8), dirmap))
            outmap = (pd.DataFrame(outmap)
                      .apply(lambda x: x.map(dir_d), axis=1).values)

        np.place(outmap, dem_mask, nodata)
        np.place(data, dem_mask, dem_nodata)

        if inplace:
            self.dir = outmap
            self.grid_props.update({'dir' : {}})
            self.grid_props['dir'].update({'bbox' : self.bbox})
            self.grid_props['dir'].update({'shape' : self.shape})
            self.grid_props['dir'].update({'cellsize' : self.cellsize})
            self.grid_props['dir'].update({'nodata' : nodata})
            self.grid_props['dir'].update({'crs' : self.crs})
        else:
            return outmap

    def catchment(self, x, y, pour_value=None, dirmap=(1, 2, 3, 4, 5, 6, 7, 8),
                  nodata=0, xytype='index', recursionlimit=15000, inplace=True):
        """
        Delineates a watershed from a given pour point (x, y).
        Returns a grid

        Parameters
        ----------
        x : int or float
            x coordinate of pour point
        y : int or float
            y coordinate of pour point
        pour_value : int or None
                     If not None, value to represent pour point in catchment
                     grid (required by some programs).
        dirmap : list or tuple
                 List of integer values representing the following
                 cardinal directions (in order):
                 [N, NE, E, SE, S, SW, W, NW]
        xytype : 'index' or 'label'
                 How to interpret parameters 'x' and 'y'.
                     'index' : x and y represent the column and row
                               indices of the pour point.
                     'label' : x and y represent geographic coordinates
                               (will be passed to self.nearest_cell).
        recursionlimit : int
                         Recursion limit--may need to be raised if
                         recursion limit is reached.
        inplace : bool
                  If True, catchment will be written to attribute 'catch'.
        """

        # if xytype is 'label', delineate catchment based on cell nearest
        # to given geographic coordinate
        if xytype == 'label':
            x, y = self.nearest_cell(x, y)

        # set recursion limit (needed for large datasets)
        sys.setrecursionlimit(recursionlimit)

        # initialize array to collect catchment cells
        self.collect = np.array([], dtype=int)

        # pad the flow direction grid with a rim of 'nodata' cells
        # easy way to prevent catchment search from going out of bounds
        try:
            self.cdir = np.pad(self.view('dir', mask=False), 1, mode='constant',
                               constant_values=np.asscalar(self.grid_props['dir']['nodata']))
        except ValueError:
            self.cdir = np.pad(self.view('dir', mask=False), 1, mode='constant')

        # get shape of padded flow direction array, then flatten
        padshape = self.cdir.shape
        self.cdir = self.cdir.ravel()

        # get the flattened index of the pour point
        pour_point = np.ravel_multi_index(np.array([y + 1, x + 1]), padshape)

        # reorder direction mapping to work with select_surround_ravel()
        dirmap = np.array(dirmap)[[4, 5, 6, 7, 0, 1, 2, 3]].tolist()

        # select cells surrounding a given cell for a flattened array
        # used by catchment_search()
        def select_surround_ravel(i):
            return np.array([i + 0 - padshape[1],
                             i + 1 - padshape[1],
                             i + 1 + 0,
                             i + 1 + padshape[1],
                             i + 0 + padshape[1],
                             i - 1 + padshape[1],
                             i - 1 + 0,
                             i - 1 - padshape[1]]).T

        # for each cell j, recursively search through grid to determine
        # if surrounding cells are in the contributing area, then add
        # flattened indices to self.collect
        def catchment_search(j):
            self.collect = np.append(self.collect, j)
            selection = select_surround_ravel(j)
            next_idx = selection[np.where(self.cdir[selection] == dirmap)]
            if next_idx.any():
                return catchment_search(next_idx)

        # call catchment search starting at the pour point
        catchment_search(pour_point)

        # initialize output array
        outcatch = np.zeros(padshape, dtype=int)

        # if nodata is not 0, replace 0 with nodata value in output array
        if nodata != 0:
            np.place(outcatch, outcatch == 0, nodata)

        # set values of output array based on 'collected' cells
        outcatch.flat[self.collect] = self.cdir[self.collect]

        # remove outer rim, delete temporary arrays
        outcatch = outcatch[1:-1, 1:-1]
        del self.cdir
        del self.collect

        # if pour point needs to be a special value, set it
        if pour_value is not None:
            outcatch[y, x] = pour_value

        # reset recursion limit
        sys.setrecursionlimit(1000)

        # if inplace is True, update attributes
        if inplace:
            self.catch = outcatch
            self.grid_props.update({'catch' : {}})
            self.grid_props['catch'].update({'bbox' : self.bbox})
            self.grid_props['catch'].update({'shape' : self.shape})
            self.grid_props['catch'].update({'cellsize' : self.cellsize})
            self.grid_props['catch'].update({'nodata' : nodata})
            self.grid_props['catch'].update({'crs' : self.crs})
        else:
            return outcatch

    def fraction(self, other, nodata=0, inplace=True):
        """
        Generates a grid representing the fractional contributing area for a
        coarse-scale flow direction grid.

        Parameters
        ----------
        other : flow_grid instance
                Another flow_grid instance containing fine-scale flow direction
                data. The ratio of self.cellsize/other.cellsize must be a
                positive integer. Grid cell boundaries must have some overlap.
                Must have attributes 'dir' and 'catch' (i.e. must have a flow
                direction grid, along with a delineated catchment).
        nodata : int or float
                 value to indicate no data in output array.

        inplace : bool (optional)
                  If True, appends fraction grid to attribute 'frac'.
        """

        # check for required attributes in self and other
        assert hasattr(self, 'dir')
        assert hasattr(other, 'dir')
        assert hasattr(other, 'catch')

        # set scale ratio
        cell_ratio = int(self.cellsize / other.cellsize)

        # create DataFrames for self and other with geographic coordinates
        # as row and column labels. entries in selfdf represent cell indices.
        selfdf = pd.DataFrame(
                np.arange(self.view('dir', mask=False).size).reshape(self.shape),
                index=np.linspace(self.bbox[1], self.bbox[3],
                                  self.shape[0], endpoint=False)[::-1],
                columns=np.linspace(self.bbox[0], self.bbox[2],
                                    self.shape[1], endpoint=False)
                )
        otherrows, othercols = self.bbox_indices(other.bbox, other.shape)

        # reindex self to other based on column labels and fill nulls with
        # nearest neighbor
        result = (selfdf.reindex(otherrows, method='nearest')
                  .reindex_axis(othercols, axis=1, method='nearest'))

        # mask cells not in catchment of 'other'
        result = result.values[np.where(other.catch != 0, True, False)]

        # count remaining indices and divide by the original number of indices
        result = ((np.bincount(result, minlength=selfdf.size)
                  .astype(float) / (cell_ratio ** 2)).reshape(selfdf.shape))

        # replace 0 with nodata value
        if nodata != 0:
            np.place(result, result == 0, nodata)

        # if inplace is True, set class attributes
        if inplace:
            self.frac = result
            self.grid_props.update({'frac' : {}})
            self.grid_props['frac'].update({'bbox' : self.bbox})
            self.grid_props['frac'].update({'shape' : self.shape})
            self.grid_props['frac'].update({'cellsize' : self.cellsize})
            self.grid_props['frac'].update({'nodata' : nodata})
            self.grid_props['frac'].update({'crs' : self.crs})
        else:
            return result

    def clip_nodata(self, data_name, precision=7, inplace=True, **kwargs):
        """
        Clip grid to bbox representing the smallest area that contains all
        non-null data for a given dataset. If inplace is True, will set
        self.bbox to the bbox generated by this method.

        Parameters
        ----------
        data_name : numpy ndarray
                    Name of attribute to base the clip on.
        inplace : bool
                  If True, update bbox to conform to clip.
        precision : int
                    Precision to use when matching geographic coordinates.
        """

        # get class attributes
        data = getattr(self, data_name)
        nodata = self.grid_props[data_name]['nodata']

        # get bbox of nonzero entries
        nz = np.nonzero(data != nodata)
        nz_ix = (nz[0].min(), nz[0].max(), nz[1].min(), nz[1].max())

        # if inplace is True, clip all grids to new bbox and set self.bbox
        if inplace:
            selfrows, selfcols = self.bbox_indices(self.grid_props[data_name]['bbox'], self.grid_props[data_name]['shape'], precision=7)
            new_bbox = (selfcols[nz_ix[2]], selfrows[nz_ix[1]],
                        selfcols[nz_ix[3]], selfrows[nz_ix[0]])
            # set self.bbox to clipped bbox
            self.set_bbox(new_bbox, **kwargs)
        else:
            # if inplace is False, return the clipped data
            return data[nz_ix[0]:nz_ix[1], nz_ix[2]:nz_ix[3]]

    def set_bbox(self, new_bbox, precision=7):
        """
        Set the bounding box of the class instance (self.bbox). If the new
        bbox is not alignable to self.cellsize, each entry  is automatically
        rounded down such that the bbox is alignable.

        Parameters
        ----------
        new_bbox : tuple
                   New bbox to use (xmin, ymin, xmax, ymax)
        precision : int
                    Precision to use when matching geographic coordinates.
        """

        # check validity of new bbox
        if ((hasattr(new_bbox, "__len__")) and (not isinstance(new_bbox, str))
                and (len(new_bbox) == 4)):
            new_bbox = tuple(new_bbox)
        else:
            raise TypeError('new_bbox must be a tuple of length 4.')

        # check if alignable; if not, round unaligned bbox entries to nearest
        new_bbox = np.asarray(new_bbox)
        err = np.abs(new_bbox) % self.cellsize
        try:
            np.testing.assert_almost_equal(err, np.zeros(len(new_bbox)))
        except AssertionError:
            off_idx = np.where(np.around(err, precision) != np.zeros(len(new_bbox)))
            direction = np.where(new_bbox > 0.0, 1, -1)
            new_bbox = new_bbox - (err * direction)
            print('Unalignable bbox provided, rounding to %s' % (new_bbox))

        # construct arrays representing old bbox coords
        selfrows, selfcols = self.bbox_indices(self.bbox, self.shape)

        # construct arrays representing coordinates of new grid
        nrows = (new_bbox[3] - new_bbox[1]) / self.cellsize
        ncols = (new_bbox[2] - new_bbox[0]) / self.cellsize
        np.testing.assert_almost_equal(nrows, round(nrows))
        np.testing.assert_almost_equal(ncols, round(ncols))
        rows = np.linspace(new_bbox[1], new_bbox[3], round(nrows), endpoint=False)
        cols = np.linspace(new_bbox[0], new_bbox[2], round(ncols), endpoint=False)

        # set class attributes
        self.bbox = tuple(new_bbox)
        self.shape = tuple([len(rows), len(cols)])
        if hasattr(self, 'catch'):
            self.catchment_mask()
        else:
            self.mask = np.ones(self.shape, dtype=np.bool)

    def set_nodata(self, data_name, new_nodata, old_nodata=None):
        """
        Change nodata value of a dataset.

        Parameters
        ----------
        data_name : string
                    Attribute name of dataset to change
        new_nodata : int or float
                     New nodata value to use
        old_nodata : int or float (optional)
                     If none provided, defaults to
                     self.grid_props[data_name]['nodata']
        """

        if old_nodata is None:
            old_nodata = self.grid_props[data_name]['nodata']
        data = getattr(self, data_name)
        np.place(data, data == old_nodata, new_nodata)
        self.grid_props[data_name]['nodata'] = new_nodata

    def catchment_mask(self, mask_source='catch'):
        """
        Masks grid cells not included in catchment. The catchment mask is saved
        to self.mask.

        Parameters
        ----------
        to_mask : string
                  name of dataset to mask
        mask_source : string (optional)
                      dataset on which mask is based (defaults to 'catch')
        """
        self.mask = (self.view(mask_source, mask=False) !=
                     self.grid_props[mask_source]['nodata'])

    def to_ascii(self, data_name=None, file_name=None, delimiter=' ', **kwargs):
        """
        Writes grid data to ascii grid files.

        Parameters
        ----------
        data_name : string or list-like (optional)
                    Attribute name(s) of datasets to write. Defaults to all
                    grid dataset names.
        file_name : string or list-like
                    Name(s) of file(s) to write to
        delimiter : string
                    Delimiter to use
        """
        if data_name is None:
            data_name = self.grid_props.keys()
        if file_name is None:
            file_name = self.grid_props.keys()

        if isinstance(data_name, str):
            data_name = [data_name]
        if isinstance(file_name, str):
            file_name = [file_name]

        for in_name, out_name in zip(data_name, file_name):
            header = """ncols         %s\nnrows         %s\nxllcorner     %s\nyllcorner     %s\ncellsize      %s\nNODATA_value  %s""" % (self.shape[1],
                                             self.shape[0],
                                             self.bbox[0],
                                             self.bbox[1],
                                             self.cellsize,
                                             self.grid_props[in_name]['nodata'])
            np.savetxt(out_name, self.view(in_name), delimiter=delimiter,
                    header=header, comments='', **kwargs)


    def _select_surround(self, i, j):
        """
        select the eight indices surrounding a given index.
        """
        return ([i - 1, i - 1, i + 0, i + 1, i + 1, i + 1, i + 0, i - 1],
                [j + 0, j + 1, j + 1, j + 1, j + 0, j - 1, j - 1, j - 1])

    def _select_edge_sur(self, edges, k):
        """
        select the five cell indices surrounding each edge cell.
        """
        i, j = edges[k]['k']
        if k == 'n':
            return ([i + 0, i + 1, i + 1, i + 1, i + 0],
                    [j + 1, j + 1, j + 0, j - 1, j - 1])
        elif k == 'e':
            return ([i - 1, i + 1, i + 1, i + 0, i - 1],
                    [j + 0, j + 0, j - 1, j - 1, j - 1])
        elif k == 's':
            return ([i - 1, i - 1, i + 0, i + 0, i - 1],
                    [j + 0, j + 1, j + 1, j - 1, j - 1])
        elif k == 'w':
            return ([i - 1, i - 1, i + 0, i + 1, i + 1],
                    [j + 0, j + 1, j + 1, j + 1, j + 0])

    def _select_surround_ravel(self, i, shape):
        """
        select the eight indices surrounding a flattened index.
        """
        return np.array([i + 0 - shape[1],
                         i + 1 - shape[1],
                         i + 1 + 0,
                         i + 1 + shape[1],
                         i + 0 + shape[1],
                         i - 1 + shape[1],
                         i - 1 + 0,
                         i - 1 - shape[1]]).T
