#!/opt/local/bin/python
"""
RVIC parameter file development driver
"""
import os
import numpy as np
import argparse
import logging
from datetime import date
from rvic.log import init_logger, log_name
from rvic.utilities import ReadConfig, MakeDirs, CopyInputs, ReadNetcdf, TarInputs
from rvic.utilities import CheckNcvars, remap, clean_file, subset, ReadDomain
from rvic.aggregate import MakeAggPairs, agg
from rvic.makeUH import rout
from rvic.share import ncGlobals
from rvic.write import write_agg_netcdf, write_param_file
from rvic.vars import Point
# try:
#     from IPython.parallel import Client
# except:
#     pass
# try:
#     import SGEpy
# except:
#     pass

# -------------------------------------------------------------------- #
# Top level driver
def main():

    UHbox, FdrData, FdrVats, DomData, Outlets, ConfigDict, dirPaths = GenUH_init()

    Outlets, Fractions, UH = GenUH_run_sp(UHbox, FdrData, FdrVats, DomData, Outlets, ConfigDict, dirPaths)

    GenUH_final(Outlets, Fractions, UH, DomData, ConfigDict, dirPaths)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Initialize the GenUH Program
def GenUH_init(configFile=None):
    """Initialize RVIC parameter """

    # ---------------------------------------------------------------- #
    # Read command Line
    if not configFile:
        configFile = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Configuration files
    ConfigDict = ReadConfig(configFile)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structure
    dirPaths = MakeDirs(ConfigDict['options']['case_dir'],
                        ['aggregated', 'remapped', 'plots',
                        'logs', 'params', 'inputs'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # copy inputs to $case_dir/inputs and update configuration
    ConfigDict = CopyInputs(configFile, dirPaths['inputs'])
    options = ConfigDict['options']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Start Logging
    log = init_logger(dirPaths['logs'], options['log_level'], options['verbose'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Pour Points files
    PourPoints = {}
    try:
        PourPoints['lons'], PourPoints['lats'] = np.genfromtxt(ConfigDict['pour_points']['file_name'],
                                                               skip_header=int(ConfigDict['pour_points']['header_lines']),
                                                               delimiter=',', unpack=True)
    except:
        log.exception('Error opening pour points file: %s' % ConfigDict['pour_points']['file_name'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read UH box file
    UHbox = {}
    try:
        UHbox['time'], UHbox['func'] = np.genfromtxt(ConfigDict['uh_box']['file_name'],
                                                     skip_header=int(ConfigDict['uh_box']['header_lines']),
                                                     delimiter=',', unpack=True)
    except:
        log.exception('Error opening uh_box file: %s' % ConfigDict['pour_points']['file_name'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read FDR file
    try:
        FdrData, FdrVats, FdrGats = ReadNetcdf(ConfigDict['routing']['file_name'])
        fdr_shape = FdrData[ConfigDict['routing']['flow_direction_var']].shape
        # Add velocity and/or diffusion grids if not present yet
        if not type(ConfigDict['routing']['velocity']) == str:
            FdrData['velocity'] = np.zeros(fdr_shape) + ConfigDict['routing']['velocity']
            ConfigDict['routing']['velocity'] = 'velocity'
            log.info('Added velocity grid to FdrData')
        if not type(ConfigDict['routing']['diffusion']) == str:
            FdrData['diffusion'] = np.zeros(fdr_shape) + ConfigDict['routing']['diffusion']
            ConfigDict['routing']['diffusion'] = 'diffusion'
            log.info('Added diffusion grid to FdrData')

        FdrData['resolution'] = np.abs(FdrData[ConfigDict['routing']['longitude_var']][1] -
                                       FdrData[ConfigDict['routing']['longitude_var']][0])

        CheckNcvars(ConfigDict['routing'], FdrData.keys())
    except:
        log.exception('Error opening FDR file')
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read domain file (if applicable)
    if options['remap']:
        DomData, DomVats, DomGats = ReadDomain(ConfigDict['domain'])
    else:
        DomData = {}
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Group pour points (if aggregate)
    if options['aggregate']:
        Outlets = MakeAggPairs(PourPoints['lons'], PourPoints['lats'],
                               DomData[ConfigDict['domain']['longitude_var']],
                               DomData[ConfigDict['domain']['latitude_var']],
                               DomData['cell_ids'], agg_type='agg')
    else:
        Outlets = {}
        for i, (lon, lat) in enumerate(zip(PourPoints['lons'], PourPoints['lats'])):
            Outlets[i] = Point(lat=lat, lon=lon)
            Outlets[i].PourPoints = [Point(lat=lat, lon=lon)]

    # ---------------------------------------------------------------- #
    return UHbox, FdrData, FdrVats, DomData, Outlets, ConfigDict, dirPaths
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# The standard single processor driver for gerating parameter files
def GenUH_run_sp(UHbox, FdrData, FdrVats, DomData, Outlets, ConfigDict, dirPaths):
    """
    Run GenUH_run on a single processor (slow)
    """
    log = logging.getLogger(log_name)

    log.info('Starting GenUH_run_sp now.')

    # ---------------------------------------------------------------- #
    # Loop over agg points
    for i, cell_id in enumerate(Outlets):

        aggData = {}
        # ------------------------------------------------------------ #
        # Loop over pour points
        for j, PourPoint in enumerate(Outlets[cell_id].PourPoints):

            # -------------------------------------------------------- #
            # Make the Unit Hydrograph Grid
            routData = rout(PourPoint, UHbox, FdrData, FdrVats,
                            ConfigDict['routing'])

            log.debug('routData: %f, %f' % (routData['UHgrid'].min(), routData['UHgrid'].max()))

            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # aggregate
            if ConfigDict['options']['aggregate']:
                if j != len(Outlets[cell_id].PourPoints)-1:
                    aggData = agg(routData, aggData, res=FdrData['resolution'])
                else:
                    aggData = agg(routData, aggData, res=FdrData['resolution'],
                                  pad=ConfigDict['options']['agg_pad'], maskandnorm=True)

                    log.debug('aggData: %f, %f' % (aggData['UHgrid'].min(), aggData['UHgrid'].max()))
            elif ConfigDict['options']['remap']:
                aggData = routData
            else:
                RmpData = routData
            # -------------------------------------------------------- #

        # ------------------------------------------------------------ #
        # write temporary file #1
        if  ConfigDict['options']['remap']:
            GlobAts = ncGlobals(title='RVIC Unit Hydrograph Grid File',
                                RvicPourPointsFile=ConfigDict['pour_points']['file_name'],
                                RvicUHFile=ConfigDict['uh_box']['file_name'],
                                RvicFdrFile=ConfigDict['routing']['file_name'],
                                RvicDomainFile=ConfigDict['domain']['file_name'])

            TempFile1 = os.path.join(dirPaths['aggregated'], 'aggUH_%i.nc' % cell_id)

            write_agg_netcdf(TempFile1, aggData, GlobAts,
                             ConfigDict['options']['netcdf_format'])
        elif ConfigDict['options']['aggregate']:
            RmpData = aggData
        else:
            pass
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Remap temporary file #1 to temporary file #2
        if ConfigDict['options']['remap']:

            TempFile2 = os.path.join(dirPaths['remapped'], 'remapUH_%i.nc' % cell_id)

            remap(ConfigDict['domain']['file_name'], TempFile1, TempFile2)

            # -------------------------------------------------------- #
            # Clean temporary file #1 (if applicable)
            if ConfigDict['options']['clean']:
                clean_file(TempFile1)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # Read temporary file #2
            RmpData, RmpVats, RmpGats = ReadNetcdf(TempFile2)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # Clean temporary file #2 (if applicable)
            if ConfigDict['options']['clean']:
                clean_file(TempFile2)
            # -------------------------------------------------------- #

        else:
            RmpData = aggData
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Add to adjust fractions Structure
        if ConfigDict['options']['remap']:
            if i == 0:
                Fractions = np.zeros(RmpData['fractions'].shape)
                UH = np.zeros(RmpData['unit_hydrographs'].shape)

            y, x = np.nonzero((RmpData['fractions'] > 0.0) * (DomData[ConfigDict['domain']['land_mask_var']] == 1))

            Outlets[cell_id].Fractions = RmpData['fractions'][y, x]
            Outlets[cell_id].UnitHydrograph = RmpData['unit_hydrographs'][:, y, x]
            Outlets[cell_id].time = np.arange(RmpData['unit_hydrographs'].shape[0])
            Outlets[cell_id].lon_source = RmpData[ConfigDict['domain']['longitude_var']][y, x]
            Outlets[cell_id].lat_source = RmpData[ConfigDict['domain']['latitude_var']][y, x]
            Outlets[cell_id].cell_id_source = DomData['cell_ids'][y, x]
            Outlets[cell_id].x_source = x
            Outlets[cell_id].y_source = y

            Fractions[y, x] += RmpData['fractions'][y, x]
            UH[:, y, x] += RmpData['unit_hydrographs'][:, y, x]
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
    return Outlets, Fractions, UH
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# The Poor Mans Parallel Driver for generating parmeter files
def GenUH_run_ppp():
    """
    Run GenUH_run in pour mans parallel (SGE only)
    """
    # ---------------------------------------------------------------- #
    # Loop over agg points
    # ---------------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # Loop over pour points
        # ------------------------------------------------------------ #
            # -------------------------------------------------------- #
            # rout(PourPoint, UH, FdrData, RoutDict)
            # -------------------------------------------------------- #
            # -------------------------------------------------------- #
            # aggregate
            # -------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # write temporary file #1
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Remap temporary file #1 to temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #1 (if applicable)
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Read temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Add to adjust fractions Structure
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #2 (if applicable)
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
    # Adjust Fractions (if applicable)
    # ---------------------------------------------------------------- #
    # Write parameter file
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def GenUH_run_ipp():
    """
    Run GenUH_run using ipython parallel
    """
    # ---------------------------------------------------------------- #
    # Loop over agg points
    # ---------------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # Loop over pour points
        # ------------------------------------------------------------ #
            # -------------------------------------------------------- #
            # RvicMakeUH()
            # -------------------------------------------------------- #
            # -------------------------------------------------------- #
            # aggregate
            # -------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # write temporary file #1
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Remap temporary file #1 to temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #1 (if applicable)
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Read temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Add to adjust fractions Structure
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #2 (if applicable)
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
    # Adjust Fractions (if applicable)
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Write parameter file
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def GenUH_final(Outlets, Fractions, UH, DomData, ConfigDict, dirPaths):
    """
    Make the RVIC Parameter File
    """

    log = logging.getLogger(log_name)

    log.info('Starting GenUH_final now.')

    options = ConfigDict['options']
    # ---------------------------------------------------------------- #
    # Determin how to adjust the fractions
    gridFracs = DomData[ConfigDict['domain']['fraction_var']]
    diffFracs = Fractions - gridFracs
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Only adjust fractions where the Fractions are gt the domain fractions
    yi, xi = np.nonzero(Fractions > gridFracs)
    ratioFracs = np.ones(diffFracs.shape)
    adjFracs = np.zeros(diffFracs.shape)
    ratioFracs[yi, xi] = gridFracs[yi, xi]/Fractions[yi, xi]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Subset
    for i, cell_id in enumerate(Outlets):

        y = Outlets[cell_id].y_source
        x = Outlets[cell_id].x_source

        offset, out_uh, full_length = subset(Outlets[cell_id].UnitHydrograph,
                                             options['subset_length'],
                                             options['subset_threshold'])
        Outlets[cell_id].Fractions *= ratioFracs[y, x]  # Adjust fracs based on ratioFracs
        adjFracs[y, x] += Outlets[cell_id].Fractions

        if i == 0:
            # -------------------------------------------------------- #
            # Source specific values
            unit_hydrographs = out_uh
            frac_sources = Outlets[cell_id].Fractions
            lon_sources = Outlets[cell_id].lon_source
            lat_sources = Outlets[cell_id].lat_source
            x_ind_source = x
            y_ind_source = y
            source_decomp_id = Outlets[cell_id].cell_id_source
            time_offset_source = offset
            source2outlet_index = np.zeros(len(Outlets[cell_id].Fractions))
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_decomp_id = np.array(cell_id)
            lon_outlet = np.array(Outlets[cell_id].lon)
            lat_outlet = np.array(Outlets[cell_id].lat)
            x_ind_outlet = np.array(Outlets[cell_id].x)
            y_ind_outlet = np.array(Outlets[cell_id].y)
            outlet_nums = np.array(i)

            # -------------------------------------------------------- #
            # get a few global values
            subset_length = options['subset_length']
            full_time_length = full_length
            unit_hydrogaph_dt = ConfigDict['routing']['output_interval']
            # -------------------------------------------------------- #

        else:
            # -------------------------------------------------------- #
            # Point specific values
            unit_hydrographs = np.append(unit_hydrographs, out_uh, axis=1)
            frac_sources = np.append(frac_sources, Outlets[cell_id].Fractions)
            lon_sources = np.append(lon_sources, Outlets[cell_id].lon_source)
            lat_sources = np.append(lat_sources, Outlets[cell_id].lat_source)
            x_ind_source = np.append(x_ind_source, x)
            y_ind_source = np.append(y_ind_source, y)
            source_decomp_id = np.append(source_decomp_id, Outlets[cell_id].cell_id_source)
            time_offset_source = np.append(time_offset_source, offset)
            source2outlet_index = np.append(source2outlet_index, np.zeros_like(offset) + i)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_decomp_id = np.append(outlet_decomp_id, cell_id)
            lon_outlet = np.append(lon_outlet, Outlets[cell_id].lon)
            lat_outlet = np.append(lat_outlet, Outlets[cell_id].lat)
            x_ind_outlet = np.append(x_ind_outlet, Outlets[cell_id].x)
            y_ind_outlet = np.append(y_ind_outlet, Outlets[cell_id].y)
            outlet_nums = np.append(outlet_nums, i)
            # -------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust Unit Hydrographs for differences in source/outlet areas
    unit_hydrographs *= frac_sources
    unit_hydrographs *= DomData[ConfigDict['domain']['area_var']][y_ind_source, x_ind_source]

    for p, ind in enumerate(source2outlet_index):
        unit_hydrographs[:, p] /= DomData[ConfigDict['domain']['area_var']][y_ind_outlet[ind], x_ind_outlet[ind]]
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Write parameter file
    today = date.today().strftime('%Y%m%d')
    ParamFile = os.path.join(dirPaths['params'],
                             '%s.rvic.prm.%s.%s.nc' % (options['caseid'],
                                                       options['gridid'], today))

    write_param_file(ParamFile, options['netcdf_format'], subset_length,
                     full_time_length, unit_hydrogaph_dt, unit_hydrographs,
                     lon_sources, lat_sources, source_decomp_id, x_ind_source,
                     y_ind_source, outlet_decomp_id, time_offset_source,
                     source2outlet_index, lon_outlet, lat_outlet, x_ind_outlet,
                     y_ind_outlet, outlet_nums,
                     ncGlobals(title='RVIC parameter file',
                               RvicPourPointsFile=os.path.split(ConfigDict['pour_points']['file_name'])[1],
                               RvicUHFile=os.path.split(ConfigDict['uh_box']['file_name'])[1],
                               RvicFdrFile=os.path.split(ConfigDict['routing']['file_name'])[1],
                               RvicDomainFile=os.path.split(ConfigDict['domain']['file_name'])[1]))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write a summary of what was done to the log file.
    log.info('Parameter file includes %i Outlets' % (len(Outlets)))
    log.info('Parameter file includes %i Source Points' % (len(Outlets)))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    InputsTar = TarInputs(dirPaths['inputs'], suffix=today)
    LogTar = TarInputs(log.filename)

    log.info('Done with RvicGenParam.')
    log.info('Location of Inputs: %s' % InputsTar)
    log.info('Location of Log: %s' % LogTar)
    log.info('Location of Parmeter File %s' % ParamFile)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def process_command_line():
    """
    Get the path to the configFile
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate RVIC parameter files.')
    parser.add_argument("configFile", type=str, help="Input configuration file")

    args = parser.parse_args()
    configFile = args.configFile

    return configFile


# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
