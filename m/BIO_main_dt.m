clear all

addpath('C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/configuration', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/setip_mex', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/tools/reading_tools', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/tools/common_tools', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/tools/netcdf_tools', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/tools/plotting_tools', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/tools/processing_tools', ...
        'C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox-master/tools/mex_tools')%, ...
%% set proper working directory
cd('C:/Users/BelzileM/Documents/Gliders/Socib/glider_toolbox/')
%% L0
configPaths = configDTPathsLocal()
configFileopts = configDTFileOptionsSeaExplorer()
configL0output = configDTOutputNetCDFL0SeaExplorer()
%   load_start = utc2posixtime(deployment_start);
%   load_final = posixtime();
%   if ~isnan(deployment_end)
%     load_final = utc2posixtime(deployment_end);
%   end
[meta, data] = loadSeaExplorerData(configPaths.ascii_path, ...              
                    configFileopts.gli_name_pattern, ...
                    configFileopts.pld_name_pattern, ...
                    'format', 'struct');
% [meta_raw, data_raw] = ...
%           loadSeaExplorerData(ascii_dir, ...
%                               file_options.gli_name_pattern, ...
%                               file_options.pld_name_pattern, ...
%                               'timegli', file_options.gli_time, ...
%                               'timepld', file_options.pld_time, ...
%                               'format', 'struct');
%         source_files = meta_raw.sources;                
ncfilename = 'C:/Users/BelzileM/Documents/Gliders/Socib/Data/SEA019/M36/nc/GLI2018_SEA019_036DM_L0.nc';
nc0 = generateOutputNetCDF(ncfilename, data, meta,'test' , ...
    configL0output.variables , configL0output.dimensions, ...
    configL0output.attributes);
% outputs.netcdf_l0 = generateOutputNetCDF( ...
%               netcdf_l0_file, data_raw, meta_raw, deployment, ...
%               netcdf_l0_options.variables, ...
%               netcdf_l0_options.dimensions, ...
%               netcdf_l0_options.attributes, ...
%               'time', {'Timestamp' 'PLD_REALTIMECLOCK'}, ...
%               'position', {'NAV_LONGITUDE' 'NAV_LATITUDE'; 'Lon' 'Lat'}, ...
%               'position_conversion', @nmea2deg, ...
%               'vertical',            {'Depth' 'SBD_PRESSURE'}, ...
%               'vertical_conversion', {[]        @(z)(z * 10)}, ...
%               'vertical_positive',   {'down'} );  
%% L1
configPrePro = configDataPreprocessingSeaExplorer();
configPro = configDataProcessingSeaExplorer()
configL1output = configDTOutputNetCDFL1()

[data_pre, meta_pre] = preprocessGliderData(data, meta, ...
                        'time_list', configPrePro.time_list, ...
                        'position_list', configPrePro.position_list, ...
                        'depth_list', configPrePro.depth_list, ...
                        'attitude_list', configPrePro.attitude_list, ...
                        'heading_list', configPrePro.heading_list, ...
                        'ctd_list', configPrePro.ctd_list, ...
                        'oxygen_list', configPrePro.oxygen_list, ...
                        'optics_list', configPrePro.optics_list, ...
                        'extra_sensor_list', configPrePro.extra_sensor_list);
% [data_preprocessed, meta_preprocessed] = ...
%         preprocessGliderData(data_raw, meta_raw, preprocessing_options);
[data_proc, meta_proc] = processGliderData(data_pre, meta_pre, ...
                            'time_filling', configPro.time_filling, ...
                            'position_filling', configPro.position_filling, ...
                            'depth_filling', configPro.depth_filling, ...
                            'attitude_filling', configPro.attitude_filling, ...
                            'heading_filling', configPro.heading_filling, ...
                            'waypoint_filling', configPro.waypoint_filling, ...
                            'pressure_filtering', configPro.pressure_filtering, ...
                            'pressure_filter_constant', configPro.pressure_filter_constant, ...
                            'depth_ctd_derivation', configPro.depth_ctd_derivation, ...
                            'profiling_list', configPro.profiling_list, ...
                            'profile_min_range', configPro.profile_min_range, ...
                            'profile_max_gap_ratio', configPro.profile_max_gap_ratio, ...
                            'sensor_lag_list', configPro.sensor_lag_list, ...
                            'thermal_lag_list', configPro.thermal_lag_list, ...
                            'salinity_list', configPro.salinity_list, ...
                            'density_list', configPro.density_list);
% [data_processed, meta_processed] = ...
%         processGliderData(data_preprocessed, meta_preprocessed, processing_options);
ncfilename = 'C:/Users/BelzileM/Documents/Gliders/Socib/Data/SEA019/M36/nc/GLI2018_SEA019_036DM_L1.nc';
nc1 = generateOutputNetCDF(ncfilename, data_proc, meta_proc,'test' , ...
    configL1output.variables , configL1output.dimensions, ...
    configL1output.attributes);
% outputs.netcdf_l1 = generateOutputNetCDF( ...
%         netcdf_l1_file, data_processed, meta_processed, deployment, ...
%         netcdf_l1_options.variables, ...
%         netcdf_l1_options.dimensions, ...
%         netcdf_l1_options.attributes);

%% L2
configGrid = configDataGridding();
configL2output = configDTOutputNetCDFL2();

[data_grid, meta_grid] = gridGliderData(data_proc, meta_proc, ...
                            'profile_list', configGrid.profile_list, ...
                            'time_list', configGrid.time_list, ...
                            'position_list', configGrid.position_list, ...
                            'depth_list', configGrid.depth_list, ...
                            'depth_step', configGrid.depth_step, ...
                            'variable_list', configGrid.variable_list);
% [data_gridded, meta_gridded] = ...
%         gridGliderData(data_processed, meta_processed, gridding_options);
ncfilename = 'C:/Users/BelzileM/Documents/Gliders/Socib/Data/SEA019/M36/nc/GLI2018_SEA019_036DM_L2.nc';
nc2 = generateOutputNetCDF(ncfilename, data_grid, meta_grid,'test' , ...
    configL2output.variables , configL2output.dimensions, ...
    configL2output.attributes);
% outputs.netcdf_l2 = generateOutputNetCDF( ...
%         netcdf_l2_file, data_gridded, meta_gridded, deployment, ...
%         netcdf_l2_options.variables, ...
%         netcdf_l2_options.dimensions, ...
%         netcdf_l2_options.attributes);
%% Figures
[figures_proc, figures_grid] = configFigures();
figure_dir='C:/Users/BelzileM/Documents/Gliders/Socib/Data/SEA019/M36/figures/';
figures = struct();
figproc_options = figures_proc;
figgrid_options = figures_grid;

% Generate processed data figures.
if ~isempty(fieldnames(data_proc)) && ~isempty(figure_dir)
    disp('Generating figures from processed data...');
    %MELANY remove complex number
    data_proc.salinity(imag(data_proc.salinity) ~= 0) = NaN;
    data_proc.chlorophyll(imag(data_proc.chlorophyll) ~= 0) = NaN;
    data_proc.density(imag(data_proc.density) ~= 0) = NaN;
    %MELANY
    try
      figures.figproc = generateGliderFigures( ...
        data_proc, figproc_options, ...
        'date', datestr(posixtime2utc(posixtime()), 'yyyy-mm-ddTHH:MM:SS+00:00'), ...
        'dirname', figure_dir);
    catch exception
      disp('Error generating processed data figures:');
      disp(getReport(exception, 'extended'));
    end
  end

% Generate gridded data figures.
  if ~isempty(fieldnames(data_grid)) && ~isempty(figure_dir)
    disp('Generating figures from gridded data...');
    try
      figures.figgrid = generateGliderFigures( ...
        data_grid, figgrid_options, ...
        'date', datestr(posixtime2utc(posixtime()), 'yyyy-mm-ddTHH:MM:SS+00:00'), ...
        'dirname', figure_dir);
    catch exception
      disp('Error generating gridded data figures:');
      disp(getReport(exception, 'extended'));
    end
  end

