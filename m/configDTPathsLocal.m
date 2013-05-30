function local_paths = configDTPathsLocal()
%CONFIGDTPATHSLOCAL  Config local paths for glider deployment delayed time data and figures.
%
%  LOCAL_PATHS = CONFIGDTPATHSLOCAL() should return a struct with the path 
%  patterns for the deployment files generated during the glider processing
%  chain in delayed time mode. It should have the following fields:
%    BINARY_PATH: path pattern of the directory containing .xbd files.
%    CACHE_PATH: path pattern of the directory for .cac files.
%    LOG_PATH: path pattern of the directory containing surface log files.
%    ASCII_PATH: path pattern of the directory for converted .dba files.
%    FIGURE_PATH: path pattern of the directory for deployment figures.
%    NETCDF_L0: path pattern of the NetCDF file for raw data
%      (data provided by the glider without any meaningful modification).
%    NETCDF_L1: path pattern of the NetCDF file for processed trajectory data
%      (properly referenced data with conversions, corrections and derivations).
%    NETCDF_L2: path pattern of the NetCDF file for processed grid data
%      (already processed data interpolated on vertical instantaneous profiles).
%    PROCESSING_LOG: path pattern of the processing log file.
%  These path patterns are converted to true paths through the function
%  STRFGLIDER.
%
%  Notes:
%    Edit this file filling in the paths to reflect your desired file layout.
%
%  Examples:
%    local_paths = configDTPathsLocal()
%
%  See also:
%    MAIN_GLIDER_DATA_PROCESSING_DT
%    STRFGLIDER
%
%  Author: Joan Pau Beltran
%  Email: joanpau.beltran@socib.cat

%  Copyright (C) 2013
%  ICTS SOCIB - Servei d'observacio i prediccio costaner de les Illes Balears.
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  error(nargchk(0, 0, nargin, 'struct'));

  local_paths.binary_path    = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/binary';
  local_paths.cache_path     = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/binary';
  local_paths.log_path       = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/log';
  local_paths.ascii_path     = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/ascii';
  local_paths.figure_path    = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/figures';
  local_paths.netcdf_l0      = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/netcdf/${GLIDER_NAME}_${DEPLOYMENT_START_DATE}_l0.nc';
  local_paths.netcdf_l1      = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/netcdf/${GLIDER_NAME}_${DEPLOYMENT_START_DATE}_l1.nc';
  local_paths.netcdf_l2      = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/netcdf/${GLIDER_NAME}_${DEPLOYMENT_START_DATE}_l2.nc';
  local_paths.processing_log = '/path/to/delayed_time/glider_data/${GLIDER_NAME}/${DEPLOYMENT_START_DATE}/processing.log';

end
