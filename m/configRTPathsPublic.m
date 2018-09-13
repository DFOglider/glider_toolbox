function public_paths = configRTPathsPublic()
%CONFIGRTPATHSPUBLIC  Configure public product and figure paths for glider deployment real time data.
%
%  Syntax:
%    PUBLIC_PATHS = CONFIGRTPATHSPUBLIC()
%
%  Description:
%    PUBLIC_PATHS = CONFIGRTPATHSPUBLIC() should return a struct
%    with the path patterns for the public copies of the deployment product
%    files generated by the glider processing chain in real time mode.
%    It should have the following fields:
%      FIGURE_DIR: path pattern of public directory for deployment figures.
%      FIGURE_URL: URL pattern pointing to public directory defined above.
%      FIGURE_INCLUDE: optional string cell array with the keys of the figures 
%        to be copied to the public location. If this fiels is not set, all
%        generated figures are copied.
%      FIGURE_EXCLUDE: optional string cell array with the keys of the figures
%        to exclude from copying to the public location.
%      FIGURE_INFO: path pattern of the public JSON file providing the list of 
%        deployment figures with their description and their URL. 
%      NETCDF_L0: path pattern of the public NetCDF file for raw data
%        (data provided by the glider without any meaningful modification).
%      NETCDF_L1: path pattern of the publict NetCDF file for processed
%        trajectory data (well referenced data with conversions, corrections,
%        and derivations).
%      NETCDF_L2: path pattern of the public NetCDF file for gridded data
%        (processed data interpolated on vertical instantaneous profiles).
%    These path patterns are converted to true paths through the function
%    STRFSTRUCT.
%
%  Notes:
%    Edit this file filling in the paths to reflect your desired file layout.
%
%  Examples:
%    public_paths = configRTPathsPublic()
%
%  See also:
%    MAIN_GLIDER_DATA_PROCESSING_RT
%    CONFIGRTPATHSLOCAL
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
 % error(nargchk(0, 0, nargin, 'struct'));
 % public_paths.figure_info = 'http://myserver/url/to/real_time/glider_data/figures/${DEPLOYMENT_ID}.json';

% pathe='c:\users\trana\matlab';
%   public_paths.netcdf_l0   = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\public\${GLIDER_NAME}_${DEPLOYMENT_START,Tyyyymmdd}_l0.nc'];
%   public_paths.netcdf_l1   = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\public\${GLIDER_NAME}_${DEPLOYMENT_START,Tyyyymmdd}_l1.nc'];
%   public_paths.netcdf_l2   = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\netcdf\${GLIDER_NAME}_${DEPLOYMENT_START,Tyyyymmdd}_l2.nc'];
%   public_paths.figure_dir  = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\public\'];
%   public_paths.figure_url  = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\public\'];
%   public_paths.figure_info = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\public\${DEPLOYMENT_ID}.json'];
%   public_paths.folder_dir  = [pathe '\glider_Project\${GLIDER_NAME}\${DEPLOYMENT_START,Tyyyymmdd}\public\'];

if ispc
    roote='w:\glider_project\';
else
    roote='/u01/rapps/glider_project/';
end
  public_paths.netcdf_l0   = fullfile(roote, 'data', 'public', '${GLIDER_NAME}_${DEPLOYMENT_START_DATE,Tyyyymmdd}_l0.nc');
  public_paths.netcdf_l1   = fullfile(roote, 'data', 'public', '${GLIDER_NAME}_${DEPLOYMENT_START_DATE,Tyyyymmdd}_l1.nc');
  public_paths.netcdf_l2   = fullfile(roote, 'data', 'public', '${GLIDER_NAME}_${DEPLOYMENT_START_DATE,Tyyyymmdd}_l2.nc');
  public_paths.figure_dir  = fullfile(roote, 'data','public');
  public_paths.figure_url  = fullfile(roote, 'data','public');
  public_paths.figure_info = fullfile(roote, 'data','public', '${DEPLOYMENT_ID}.json');
  public_paths.folder_dir  = fullfile(roote, 'data','public');

end
