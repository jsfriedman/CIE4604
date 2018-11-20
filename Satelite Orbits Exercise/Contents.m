% TLE (Satellite NORAD Two Line Elements) Toolbox
% Version 1.2-l (September 13, 2017) - Limited edition for CIE4604
%
% Get and read NORAD Two Line Elements in TLE structure array
%   tleget      - Retrieve NORAD Two Line Elements from www.celestrak.com
%   tleread     - Read NORAD Two Line Elements from file.
%
% Plotting of satellite positions, elevation, azimuth, visibility, rise/set, ...
%   tleplot1    - Plot satellite position and velocity from NORAD Two Line Elements.
%
% Find satellites and select dates 
%   tlefind     - Find named satellites in the NORAD Two Line Elements.
%   tledatenum  - Compute Matlab datenumbers from a date range.
%
% Compute satellite positions and orbit propagation
%   tle2vec1    - Satellite position and velocity from NORAD Two Line Elements.
%   tle2orb     - Compute orbital elements from NORAD Two Line Elements
%
% Examples:
%   tle=tleread('resource-10-oct-2017.tle');
%   tlefind(tle,'SENTINEL');
%   tleplot1(tle,{'2017-10-10 0:00', 24*60 ,1},'SENTINEL-1A',[ 52 4.8  0 ]); 
%
% This toolbox requires functions from the CRSUTIL toolbox.
%
% See also CRSUTIL.
%
% (c) Hans van der Marel, Delft University of Technology, 2012-2017.
