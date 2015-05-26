import numpy as np
import pandas as pd
import sys
import ast

class flow_grid(object):        
    """
    Container class for holding and manipulating gridded VIC routing data.
    Can be instantiated with optional keyword arguments. These keyword
    arguments will be used to add a dataset (dem, flowdir, accumulation)
    to the flow_grid instance.

    Parameters
    ----------
    data : string or numpy ndarray (optional)
           Data to be read. Can either be a file name or an array.
           If data is from a file, 'input_type' should be set to the
           appropriate value ('ascii' or 'raster').
    data_type : 'dem', 'dir', 'acc' (optional)
                 How to interpret the input data:
                     'dem' : digital elevation data
                     'dir' : flow direction data
                     'acc' : flow accumulation (upstream area) data
                
    input_type : 'raster', 'ascii' or 'array' (optional)
                 Type of input data.
    band : int (optional)
           For raster data, the band number to read.
    nodata : int or float (optional)
             Value indicating no data.

    Attributes (Optional)
    ---------------------
    dem : digital elevation grid
    dir : flow direction grid
    acc : flow accumulation grid
    catch : Catchment delineated from 'dir' and a given pour point
    frac : fractional contributing area grid
    bbox : The geographical bounding box of the gridded dataset
           (xmin, ymin, xmax, ymax)
    shape : The shape of the gridded data (nrows, ncolumns)
    cellsize : The length/width of each grid cell (assumed to be square).
    nodata : The value to use for gridcells with no data.

    Methods
    -------
    read_input : add a gridded dataset (dem, flowdir, accumulation) 
                 to flow_grid instance.
    nearest_cell : Returns the index (column, row) of the cell closest
                   to a given geographical coordinate (x, y).
    flowdir : Generate a flow direction grid from a given digital elevation
              dataset (dem).
    catchment : Delineate the watershed for a given pour point (x, y)
                or (column, row).
    fraction : Generate the fractional contributing area for a coarse
               scale flow direction grid based on a fine-scale flow
               direction grid.
    """

    def __init__(self, conform_data=True, **kwargs):
        self.conform_data = conform_data
        self.nodata = {}
        self.gridlist = []
        if 'data' in kwargs:
            self.read_input(**kwargs)
        else:
            pass

    def read_input(self, data, data_type='dir', input_type='ascii', band=1,
                   nodata=0, bbox=None, crs=None, **kwargs):
        """
        Reads data into a named attribute of flow_grid
        (name of attribute determined by 'data_type').

        Parameters
        ----------
        data : File name (string) or numpy ndarray
               If data is from a file, 'input_type' should be set to the
               appropriate value ('ascii' or 'raster').
        data_type : 'dem', 'dir', 'acc' or other string
                     Name of dataset. Will determine the name of the attribute
                     representing the gridded data. Default values are used
                     internally by some class methods:
                         'dem' : digital elevation data
                         'dir' : flow direction data
                         'acc' : flow accumulation (upstream area) data
                    
        input_type : 'raster', 'ascii' or 'array'
                     Type of input data.
        band : int
               For raster data, the band number to read.
        nodata : int or float
                 Value indicating no data.
        bbox : tuple or list
               Bounding box, if none provided.

        """

        # read ascii file
        if input_type == 'ascii':
            with open(data) as header:
                ncols = int(header.readline().split()[1])
                nrows = int(header.readline().split()[1])
                xll = ast.literal_eval(header.readline().split()[1])
                yll = ast.literal_eval(header.readline().split()[1])
                cellsize = ast.literal_eval(header.readline().split()[1])
                nodata = ast.literal_eval(header.readline().split()[1])
                shape = (nrows, ncols)
                bbox = (xll, yll, xll + ncols*cellsize, yll + nrows*cellsize)
            data = np.loadtxt(data, skiprows=6, **kwargs)
            nodata = data.dtype.type(nodata)

        # read raster file
        if input_type == 'raster':
            import rasterio
            f = rasterio.open(data)
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

        # read numpy array
        if input_type == 'array':
            shape = data.shape

        # if there are no datasets, initialize bbox, shape,
        # cellsize and crs based on incoming data
        if len(self.gridlist) < 1:
            self.bbox = bbox
            self.shape = shape
            self.cellsize = cellsize
            self.crs = crs
        # if there are existing datasets, conform incoming
        # data to bbox
        else:
            if self.conform_data == True:
                np.testing.assert_almost_equal(cellsize, self.cellsize)
                try:
                    np.testing.assert_almost_equal(bbox, self.bbox)
                except AssertionError:
                    data = self.conform(data, bbox, fill=nodata)
                    shape = data.shape
                assert(shape == self.shape)

        self.shape_min = np.min_scalar_type(max(self.shape))
        self.size_min = np.min_scalar_type(data.size)

        # assign new data to attribute; record nodata value
        self.gridlist.append(data_type)
        self.nodata.update({data_type : nodata})
        setattr(self, data_type, data)

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
            np.linspace(self.bbox[0] + self.cellsize/2.0,
                        self.bbox[2] + self.cellsize/2.0,
                        self.shape[1], endpoint=False),
            np.linspace(self.bbox[1] + self.cellsize/2.0,
                        self.bbox[3] + self.cellsize/2.0,
                        self.shape[0], endpoint=False)[::-1])

        # select nearest cell based on euclidian distance
        nearest = np.unravel_index(np.argmin(np.sqrt((
                                   coords[0] - x)**2 + (coords[1] - y)**2)),
                                   self.shape)
        return nearest[1], nearest[0]

    def flowdir(self, data=None, include_edges=True, dirmap=(1,2,3,4,5,6,7,8),
                nodata=0, flat=-1, inplace=True):
        """
        Generates a flow direction grid from a DEM grid.

        Parameters
        ----------
        data : numpy ndarray
               Array representing DEM grid
        include_edges : bool
                        Whether to include outer rim of grid.
        dirmap : list or tuple
                 List of integer values representing the following
                 cardinal and intercardinal directions (in order):
                 [N, NE, E, SE, S, SW, W, NW]
        """

        # if data not provided, use self.dem
        if data is None:
            if hasattr(self, 'dem'):
                data = self.dem

        # generate grid of indices
        indices = np.indices(self.shape, dtype=self.shape_min)

        # handle nodata values in dem
        dem_nodata = self.nodata['dem']
        dem_mask = (data == dem_nodata)        
        np.place(data, dem_mask, np.iinfo(data.dtype.type).max)

        # initialize indices of corners
        corner = {
        'nw' : {'k' : tuple(indices[:,0,0]),
                'v' : [[0,1,1], [1,1,0]],
                'pad': np.array([3,4,5])},
        'ne' : {'k' : tuple(indices[:,0,-1]),
                'v' : [[1,1,0], [-1,-2,-2]],
                'pad': np.array([5,6,7])},
        'sw' : {'k' : tuple(indices[:,-1,0]),
                'v' : [[-2,-2,-1], [0,1,1]],
                'pad': np.array([1,2,3])},
        'se' : {'k' : tuple(indices[:,-1,-1]),
                'v' : [[-1,-2,-2], [-2,-2,-1]],
                'pad': np.array([7,8,1])}
        }
    
        # initialize indices of edges
        edge = {
        'n' : {'k' : tuple(indices[:,0,1:-1]),
               'pad' : np.array([3,4,5,6,7])},
        'w' : {'k' : tuple(indices[:,1:-1,0]),
               'pad' : np.array([1,2,3,4,5])},
        'e' : {'k' : tuple(indices[:,1:-1,-1]),
               'pad' : np.array([1,5,6,7,8])},
        's' : {'k' : tuple(indices[:,-1,1:-1]),
               'pad' : np.array([1,2,3,7,8])}
        }
    
        # initialize indices of body (all cells except edges and corners)
        body = indices[:, 1:-1, 1:-1]
    
        # initialize output array
        outmap = np.full(self.shape, nodata, dtype=np.int8)
    
        # select the eight cells surrounding a cell
        def select_surround(i, j):
            return ([i-1, i-1, i+0, i+1, i+1, i+1, i+0, i-1],
                   [j+0, j+1, j+1, j+1, j+0, j-1, j-1, j-1])
    
        # select the five cells surrounding an edge cell
        def select_edge_sur(k):
            i,j = edge[k]['k']
            if k == 'n':
                return [i+0, i+1, i+1, i+1, i+0], [j+1, j+1, j+0, j-1, j-1]
            elif k =='e':
                return [i-1, i+1, i+1, i+0, i-1], [j+0, j+0, j-1, j-1, j-1]
            elif k =='s':
                return [i-1, i-1, i+0, i+0, i-1], [j+0, j+1, j+1, j-1, j-1]
            elif k == 'w':
                return [i-1, i-1, i+0, i+1, i+1], [j+0, j+1, j+1, j+1, j+0]
     
        # for each entry in "body" determine flow direction based
        # on steepest neighboring slope
        for i, j in np.nditer(tuple(body), flags=['external_loop']):
            dat = data[i,j]
            sur = data[select_surround(i,j)]
            a = ((dat - sur) > 0).any(axis=0)
            b = np.argmax((dat - sur), axis=0) + 1
            c = flat
            outmap[i,j] = np.where(a,b,c)

        # determine flow direction for edges and corners, if desired
        if include_edges == True:

            # fill corners
            for i in corner.keys():
                dat = data[corner[i]['k']]
                sur = data[corner[i]['v']]
                if ((dat - sur) > 0).any():
                    outmap[corner[i]['k']] = corner[i]['pad'][np.argmax(dat - sur)]
                else:
                    outmap[corner[i]['k']] = flat

            # fill edges
            for x in edge.keys():
                dat = data[edge[x]['k']]
                sur = data[select_edge_sur(x)]
                a = ((dat - sur) > 0).any(axis=0)
                b = edge[x]['pad'][np.argmax((dat - sur), axis=0)]
                c = flat
                outmap[edge[x]['k']] = np.where(a,b,c)
    
        # If direction numbering isn't default, convert values of output array.
        if dirmap != (1,2,3,4,5,6,7,8):
            dir_d = dict(zip((1,2,3,4,5,6,7,8), dirmap))
            outmap = pd.DataFrame(outmap).apply(lambda x: x.map(dir_d), axis=1).values

        np.place(outmap, dem_mask, nodata)
        np.place(data, dem_mask, dem_nodata)        
        
        if inplace == True:
            self.dir = outmap
        else:
            return outmap

    def catchment(self, x, y, pour_value=None, dirmap=(1,2,3,4,5,6,7,8),
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
            self.cdir = np.pad(self.dir, 1, mode='constant',
                               constant_values=np.asscalar(self.nodata['dir']))
        except ValueError:
            self.cdir = np.pad(self.dir, 1, mode='constant')

        # get shape of padded flow direction array, then flatten
        padshape = self.cdir.shape
        self.cdir = self.cdir.ravel()
        
        # get the flattened index of the pour point
        pour_point = np.ravel_multi_index(np.array([y+1, x+1]), padshape)

        # reorder direction mapping to work with select_surround_ravel()
        dirmap = np.array(dirmap)[[4,5,6,7,0,1,2,3]].tolist()

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
            outcatch[y,x] = pour_value 

        # if inplace is True, update attributes
        if inplace == True:
            self.catch = outcatch
            self.nodata.update({'catch' : nodata})
            self.gridlist.append('catch')
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
                
        inplace : bool (optional)
                  If True, appends fraction grid to attribute 'frac'.
        """

        # check for required attributes in self and other
        assert hasattr(self, 'dir')
        assert hasattr(other, 'dir')
        assert hasattr(other, 'catch')

        # set scale ratio
        cell_ratio = int(self.cellsize/other.cellsize)

        # create DataFrames for self and other with geographic coordinates
        # as row and column labels. entries in selfdf represent cell indices.
        selfdf = pd.DataFrame(
                np.arange(self.dir.size).reshape(self.shape),
                index=np.linspace(self.bbox[1], self.bbox[3],
                                  self.shape[0], endpoint=False)[::-1],
                columns=np.linspace(self.bbox[0], self.bbox[2],
                                    self.shape[1], endpoint=False)
                )
        otherdf = pd.DataFrame(
                other.dir,
                index=np.linspace(other.bbox[1], other.bbox[3],
                                  other.shape[0], endpoint=False)[::-1],
                columns=np.linspace(other.bbox[0], other.bbox[2],
                                    other.shape[1], endpoint=False)
                )

        # reindex self to other based on column labels and fill nulls with
        # nearest neighbor
        result = (selfdf.reindex(otherdf.index, method='nearest')
                  .reindex_axis(otherdf.columns, axis=1, method='nearest'))

        # mask cells not in catchment of 'other'
        result = result.values[np.where(other.catch != 0, True, False)]

        # count remaining indices and divide by the original number of indices
        result = ((np.bincount(result, minlength=selfdf.size)
                  .astype(float)/(cell_ratio**2)).reshape(selfdf.shape))

        # replace 0 with nodata value
        if nodata != 0:
            np.place(result, result == 0, nodata)

        # if inplace is True, set class attributes
        if inplace == True:
            self.frac = result
            self.nodata.update({'frac' : nodata})
            self.gridlist.append('frac')
        else:
            return result

    def conform(self, data, bbox, precision=7, fillna=True, fill=0):
        """
        Conform data to existing bbox.

        Parameters
        ----------
        data : numpy ndarray
               Data to be conformed.
        bbox : tuple
               bbox of new data
        precision : int
                    Precision to use when matching geographic coordinates.
        fillna : bool
                 Whether to fill nulls.
        fill : int or float
               Fill value to use.
        """

        # create arrays representing coordinates of existing grids
        selfrows = np.around(np.linspace(self.bbox[1], self.bbox[3],
                   self.shape[0], endpoint=False)[::-1], precision)
        selfcols = np.around(np.linspace(self.bbox[0], self.bbox[2],
                   self.shape[1], endpoint=False), precision)

        # create arrays representing coordinates of new grid
        rows = np.around(np.linspace(bbox[1], bbox[3], data.shape[0],
                                     endpoint=False)[::-1], precision)
        cols = np.around(np.linspace(bbox[0], bbox[2], data.shape[1],
                                     endpoint=False), precision)

        # reindex new grid to existing grids
        data = pd.DataFrame(data, index=rows, columns=cols)
        data = data.reindex(selfrows).reindex_axis(selfcols, axis=1)

        # if there is area with no overlap, fill nulls
        if fillna == True:
            return data.fillna(fill).values
        else:
            return data.values

    def clip_nodata(self, data_name, inplace=True, precision=7, **kwargs):
        """
        Clip grid to bbox representing the smallest area that contains all
        non-null data. If inplace is True, will also set self.bbox to this
        value and clip all other grids to the same extent.

        Parameters
        ----------
        data_name : numpy ndarray
                    Name of attribute to base the clip on.
        precision : int
                    Precision to use when matching geographic coordinates.
        """

        # get class attributes
        data = getattr(self, data_name)
        nodata = self.nodata[data_name]

        # get bbox of nonzero entries
        nz = np.nonzero(data != nodata)
        nz_ix = (nz[0].min(), nz[0].max(), nz[1].min(), nz[1].max())

        # if inplace is True, clip all grids to new bbox and set self.bbox
        if inplace == True:
            selfrows = np.around(np.linspace(self.bbox[1], self.bbox[3],
                       self.shape[0], endpoint=False)[::-1],precision)
            selfcols = np.around(np.linspace(self.bbox[0], self.bbox[2],
                       self.shape[1], endpoint=False), precision)
            new_bbox = (selfcols[nz_ix[2]], selfrows[nz_ix[1]],
                        selfcols[nz_ix[3]], selfrows[nz_ix[0]])
            # set self.bbox to clipped bbox
            self.set_bbox(new_bbox, **kwargs)
        else:
            # if inplace is False, return the clipped data
            return data[nz_ix[0]:nz_ix[1], nz_ix[2]:nz_ix[3]]

    def set_bbox(self, new_bbox, precision=7): 
        """
        Set the bounding box of the class instance (self.bbox) and clip
        all grids to the new bbox.

        Parameters
        ----------
        new_bbox : tuple
                   New bbox to use (xmin, ymin, xmax, ymax)
        precision : int
                    Precision to use when matching geographic coordinates.
        fillna : bool
                 Whether to fill nulls.
        fill : int or float
               Fill value to use.
        """
        # round new bbox to proper precision
        new_bbox = np.around(new_bbox, precision)

        # construct arrays representing coordinates of existing grids
        selfrows = np.around(np.linspace(self.bbox[1], self.bbox[3],
                   self.shape[0], endpoint=False)[::-1], precision)
        selfcols = np.around(np.linspace(self.bbox[0], self.bbox[2],
                   self.shape[1], endpoint=False), precision)

        # construct arrays representing coordinates of new grid
        rows = (pd.Series(selfrows, index=selfrows)
                .loc[new_bbox[3]:new_bbox[1]].values)
        cols = (pd.Series(selfcols, index=selfcols)
                .loc[new_bbox[0]:new_bbox[2]].values)

        # clip existing grids to new bbox
        for i in self.gridlist:
            data = pd.DataFrame(getattr(self, i),
                                index=selfrows, columns=selfcols)
            data = (data.reindex(rows).reindex_axis(cols, axis=1)
                        .fillna(self.nodata[i]).values)
            setattr(self, i, data)
        
        # set class attributes
        self.bbox = tuple(new_bbox)
        self.shape = tuple([len(rows), len(cols)])

    def set_nodata(self, data_name, new_nodata, old_nodata=None):
        """
        Change nodata value of dataset.

        Parameters
        ----------
        data_name : string
                    Attribute name of dataset to change
        new_nodata : int or float
                     New nodata value to use
        old_nodata : int or float (optional)
                     If none provided, defaults to self.nodata[data_name]
        """

        if old_nodata is None:
            old_nodata = self.nodata[data_name]
        data = getattr(self, data_name)
        np.place(data, data==old_nodata, new_nodata)
        self.nodata.update({data_name : new_nodata})

    def catchment_mask(self, to_mask, mask_source='catch'):
        """
        Masks grid cells not included in catchment.

        Parameters
        ----------
        to_mask : string
                  name of dataset to mask
        mask_source : string (optional)
                      dataset on which mask is based (defaults to 'catch')
        """

        mask = (getattr(self, mask_source) == self.nodata[mask_source])
        if isinstance(to_mask, str):
            np.place(getattr(self, to_mask), mask, self.nodata[to_mask])
        elif isinstance(to_mask, (list, tuple, np.ndarray)):
            for i in to_mask:
                np.place(getattr(self, i), mask, self.nodata[i])

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
            data_name = self.gridlist
        if file_name is None:
            file_name = self.gridlist

        if isinstance(data_name, str):
            data_name = [data_name]
        if isinstance(file_name, str):
            file_name = [file_name]

        for i, j in zip(data_name, file_name):
            header = """ncols         %s\n
                        nrows         %s\n
                        xllcorner     %s\n
                        yllcorner     %s\n
                        cellsize      %s\n
                        NODATA_value  %s\n""" % (self.shape[1],
                                               self.shape[0],
                                               self.bbox[0],
                                               self.bbox[1],
                                               self.cellsize,
                                               self.nodata[i])
            np.savetxt(j, getattr(self, i),
                       delimiter=delimiter, header=header, **kwargs)
