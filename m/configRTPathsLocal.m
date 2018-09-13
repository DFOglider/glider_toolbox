function local_paths = configRTPathsLocal()
%CONFIGRTPATHSLOCAL  Config local paths for glider deployment real time data and figures.
%
%  Syntax:
%    LOCAL_PATHS = CONFIGRTPATHSLOCAL()
%
%  Description:
%    LOCAL_PATHS = CONFIGDTPATHSLOCAL() should return a struct with the path 
%    patterns for the deployment files involved in the glider processing chain
%    in real time mode. It should have the following fields:
%      BINARY_PATH: path pattern of directory of Slocum binary files (.Xbd).
%      CACHE_PATH: path pattern of directory of Slocum cache files (.cac).
%      LOG_PATH: path pattern of directory of Slocum surface log files (.log).
%      ASCII_PATH: path pattern of directory of converted Slocum text files
%        (.dba) or Seglider text files (.log and .eng).
%      FIGURE_PATH: path pattern of directory for deployment figures.
%      NETCDF_L0: path pattern of NetCDF file for raw data
%        (data provided by the glider without any meaningful modification).
%      NETCDF_L1: path pattern of NetCDF file for processed trajectory data
%        (well referenced data with conversions, corrections and derivations).
%      NETCDF_L2: path pattern of NetCDF file for processed grid data
%        (processed data interpolated on vertical instantaneous profiles).
%      PROCESSING_LOG: path pattern of processing log file.
%    These path patterns are converted to true paths through the function
%    STRFSTRUCT.
%
%  Notes:
%    Edit this file filling in the paths to reflect your desired file layout.
%
%  Examples:
%    local_paths = configRTPathsLocal()
%
%  See also:
%    MAIN_GLIDER_DATA_PROCESSING_RT
%    STRFSTRUCT
%
%  Authors:
%    Joan Pau Beltran  <joanpau.beltran@socib.cat>

%  Copyright (C) 2013-2016
%  ICTS SOCIB - Servei d'observacio i prediccio costaner de les Illes Balears
%  <http://www.socib.es>
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
narginchk(0, 0);
%  error(nargchk(0, 0, nargin, 'struct'));
if ispc
    roote='w:\glider_project\';
else
    roote='/u01/rapps/glider_project/';
end

%  % pathe='c:\users\trana\matlab';
%   local_paths.binary_path    = [roote '\data\${GLIDER_TYPE}\input\'];
%   local_paths.cache_path     = [roote '\data\${GLIDER_TYPE}\input\cache\'];
%   local_paths.log_path       = [roote '\data\log\'];
%   local_paths.ascii_path     = [roote '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\ascii\'];
%   local_paths.figure_path    = [roote '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\output\'];
%   local_paths.netcdf_l0      = [roote '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\output\${GLIDER_NAME}_${DEPLOYMENT_START,Tyyyymddm}_l0.nc'];
%   local_paths.netcdf_l1      = [roote '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\output\${GLIDER_NAME}_${DEPLOYMENT_START,Tyyyymmdd}_l1.nc'];
%   local_paths.netcdf_l2      = [roote '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\output\${GLIDER_NAME}_${DEPLOYMENT_START,Tyyyymmdd}_l2.nc'];
%   local_paths.processing_log = [roote '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\output\processing.log'];
  local_paths.binary_path    = fullfile(roote, 'data', '${GLIDER_MODEL}', 'input');
  local_paths.cache_path     = fullfile(roote, 'data', '${GLIDER_MODEL}', 'input', 'cache');
  local_paths.log_path       = fullfile(roote ,'data', 'log');
  local_paths.ascii_path     = fullfile(roote, 'data', '${GLIDER_MODEL}', 'input');
  local_paths.figure_path    = fullfile(roote ,'data', '${GLIDER_MODEL}','output');  
%     local_paths.netcdf_l0      = fullfile(roote, 'data', '${GLIDER_MODEL}', 'output', '${GLIDER_NAME}_${DEPLOYMENT_ID}_l0.nc');
%     local_paths.netcdf_l1      = fullfile(roote, 'data', '${GLIDER_MODEL}', 'output', '${GLIDER_NAME}_${DEPLOYMENT_ID}_l1.nc');
%     local_paths.netcdf_l2      = fullfile(roote, 'data', '${GLIDER_MODEL}', 'output', '${GLIDER_NAME}_${DEPLOYMENT_ID}_l2.nc');
    local_paths.netcdf_l0      = fullfile(roote, 'data', '${GLIDER_MODEL}', 'output', '${GLIDER_NAME}_${DEPLOYMENT_START_DATE,Tyyyymmdd}_l0.nc');
    local_paths.netcdf_l1      = fullfile(roote, 'data', '${GLIDER_MODEL}', 'output', '${GLIDER_NAME}_${DEPLOYMENT_START_DATE,Tyyyymmdd}_l1.nc');
    local_paths.netcdf_l2      = fullfile(roote, 'data', '${GLIDER_MODEL}', 'output', '${GLIDER_NAME}_${DEPLOYMENT_START_DATE,Tyyyymmdd}_l2.nc');
  local_paths.processing_log = fullfile(roote, 'data', 'log','${GLIDER_NAME}_${DEPLOYMENT_ID}_processing.log');
end
