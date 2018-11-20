function tle=tleread(fname,varargin)
%TLEREAD   Read NORAD Two Line Elements from file.
%   TLE=TLEREAD(FNAME) reads NORAD two line elements from file and output
%   the TLE elements in the structure array TLE.
%
%   The mean orbital elements contained in TLE are Earth-centered-inertial (ECI)
%   coordinates with respect to the true equator of date and the mean equinox 
%   of date. They do not include the effect of nutation. The TLE are only 
%   compatible with the SGP4 and SDP4 orbit propagators of NORAD and should
%   not be used with other orbit propagators.
%
%   The TLE files can have an optional line with a twenty-four character name 
%   before the traditional Two Line Element format (Three-Line Elemement Set). 
%
%   Files with TLE's can be obtained from www.celestrak.com. You may use
%   the function TLEGET to do this.
%
%   Example:
%     tle=tleread('gps-ops.tle')
%
%   See also TLEGET, TLE2VEC and TLE2AZEL.

%   (c) Hans van der Marel, Delft Universtiy of Technology, 2012-2013.
%
%   Created:    30 Dec 2012 by Hans van der Marel
%   Modified:    4 Jan 2013 by Hans van der Marel
%                  - added support for TLE files without line0
%                2 Nov 2013 by Hans van der Marel
%                  - changed lay-out of verbose print 
%                3 Aug 2015 by Hans van der Marel
%                  - check for existence of fname
%               25 Aug 2015 by Hans van der Marel
%                  - support for option pairs (verbose)
%               13 Sep 2017 by Hans van der Marel
%                  - change of tle structure format

% Constants (WGS84)

mu = 398600.5;           %  Earth gravitational parameter (WGS84) [ km^3 / s^2 ]
Re = 6378.137;           %  Earth radius (WGS84) [ km ]
d2r = pi/180;            %  Degrees to radians

% Check the options

opt.verbose=1;

for k=1:2:length(varargin)
   if isfield(opt,varargin{k})
     opt.(varargin{k})=varargin{k+1};
   else
     error(['Invalid option ' varargin{k}])
   end
end

% Open the TLE file

if ~exist(fname,'file')
   error(['Filename '  fname ' with TLE elements not found.']);
end

fid = fopen(fname);

% Data for each satellite consists of three lines in the following format:
%
%          1         2         3         4         5         6         7
% 1234567890123456789012345678901234567890123456789012345678901234567890
%
% AAAAAAAAAAAAAAAAAAAAAAAA
% 1 NNNNNU NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
% 2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
%
% Line 0 is a twenty-four character name (to be consistent with the name 
% length in the NORAD SATCAT). Line 0 is optional and may be preceeded with
% a zero.  
%
% Lines 1 and 2 are the standard Two-Line Orbital Element Set Format 
% identical to that used by NORAD and NASA. The format description is:
%
% Line 1
%
% Column 	Description
% 01        Line Number of Element Data
% 03-07 	Satellite Number
% 08        Classification (U=Unclassified)
% 10-11 	International Designator (Last two digits of launch year)
% 12-14 	International Designator (Launch number of the year)
% 15-17 	International Designator (Piece of the launch)
% 19-20 	Epoch Year (Last two digits of year)
% 21-32 	Epoch (Day of the year and fractional portion of the day)
% 34-43 	First Time Derivative of the Mean Motion
% 45-52 	Second Time Derivative of Mean Motion (decimal point assumed)
% 54-61 	BSTAR drag term (decimal point assumed)
% 63        Ephemeris type
% 65-68 	Element number
% 69        Checksum (Modulo 10)
%
% Line 2
%
% Column 	Description
% 01        Line Number of Element Data
% 03-07 	Satellite Number
% 09-16 	Inclination [Degrees]
% 18-25 	Right Ascension of the Ascending Node [Degrees]
% 27-33 	Eccentricity (decimal point assumed)
% 35-42 	Argument of Perigee [Degrees]
% 44-51 	Mean Anomaly [Degrees]
% 53-63 	Mean Motion [Revs per day]
% 64-68 	Revolution number at epoch [Revs]
% 69        Checksum (Modulo 10)
% 
% All other columns are blank or fixed.
%
% The mean orbital elements in TLE are Earth-centered-inertial (ECI)
% coordinates with respect to the true equator of date and the mean equinox 
% of date. They do not include the effect of nutation.


% Read TLE elements

ntle=0;
%if opt.verbose, fprintf('\n a [km]      ecc      inc [deg] RAAN [deg] argp [deg]   E [deg]   Period     Reference_Epoch     Satellite\n');, end
if opt.verbose, fprintf('\nSatellite                  Reference_Epoch    a [km]      ecc      inc [deg] RAAN [deg] argp [deg]   E [deg]   Period\n');, end
while ~feof(fid)

   % Decode the line with the satellite name (optional) and read next two lines
   line0=fgetl(fid);
   if strcmp(line0(1:1),'1')
      threeline=0;
   else 
      threeline=1;
   end  
   if threeline 
      if strcmp(line0(1:1),'0')
         satname=sscanf(line0(3:end),'%24c%*s',1);
      else
         satname=sscanf(line0,'%24c%*s',1);
      end
      line1=fgetl(fid);
   else
      [~,satname] = fileparts(fname);
      line1=line0;
   end
   if  sscanf(line1(1:1),'%d') ~= 1  
      error(['Incorrect line number TLE ' satname ', 1 expected: ' line1 ])
   end
   line2=fgetl(fid);
   if ( sscanf(line2(1:1),'%d') ~= 2 )
      error(['Incorrect line number TLE ' satname ', 2 expected: ' line2 ])
   end 

   % Decode the first line with TLE data
   satid=sscanf(line1(3:7),'%d');        % Satellite Number
   classification=line1(8:8);            % Classification (U=Unclassified)
   intldesg = line1(10:17);              % International Designator 
   epochyr = sscanf(line1(19:20),'%d');  % Epoch Year (Last two digits of year)
   if (epochyr < 57)
       epochyr=epochyr + 2000;
   else
       epochyr=epochyr + 1900;
   end
   epochdays =sscanf(line1(21:32),'%f'); % Epoch (Day of the year and fractional portion of the day)
   ndot = sscanf(line1(34:43),'%f');     % First Time Derivative of the Mean Motion [rev/day^2]
   nddot = sscanf(line1(45:50),'%d');    % Second Time Derivative of Mean Motion  [rev/day^3] (decimal point assumed)
   nexp = sscanf(line1(51:52),'%d');
   nddot = nddot*1e-5*10.0^nexp;
   bstar = sscanf(line1(54:59),'%d');    % BSTAR drag term (decimal point assumed)
   ibexp = sscanf(line1(60:61),'%d');
   bstar = bstar*1e-5*10.0^ibexp;
   ephtype =line1(63:63);                % Ephemeris type
   elnum = sscanf(line1(65:68),'%d');    % Element number   

   % Decode the second line with TLE data
   if ( sscanf(line2(3:7),'%d') ~= satid )
      error(['Satellite id on 2nd TLE does not match first ' satid ', line 2 ' line2 ])
   end
   inc0 = sscanf(line2(8:16),'%f');         % Inclination [deg]
   raan0 = sscanf(line2(17:25),'%f');       % Right Ascension of the Ascending Node [deg]
   ecc0 = sscanf(line2(27:33),'%d')*1e-7;   % Eccentricity (decimal point assumed)
   argp0 = sscanf(line2(34:42),'%f');       % Argument of periapsis [deg]
   m0 = sscanf(line2(43:51),'%f');          % Mean anomaly [deg]
   n0 = sscanf(line2(52:63),'%f');          % Mean motion [Revs per day]
   revnum = sscanf(line2(64:68),'%d');      % Revolution number [Revs]
   
   % Complete orbital elements
   t0=datenum(epochyr,1,floor(epochdays))+epochdays-floor(epochdays);
   a0 = (mu/(n0*2*pi/(24*3600))^2)^(1/3);   % Semi-major axis [km]    
   e0 = eanomaly(m0,ecc0);                  % Eccentric anomaly
   OE = [a0 ecc0 inc0 raan0 argp0 e0];
   if opt.verbose
      ihour=floor(24./n0);
      imin=floor(24*60./n0-ihour*60);
      isec=round(24*3600./n0-ihour*3600-imin*60);
      if isec==60
          isec=0;
          imin=imin+1;
      end
      tt=sprintf('%02d:%02d:%02d',ihour,imin,isec);
      %fprintf('%8.2f  %9.7f  %9.4f  %9.4f  %9.4f  %9.4f  %s  %s  %s\n', OE,tt,datestr(t0,0),satname);
      fprintf('%24s%s %8.2f  %9.7f  %9.4f  %9.4f  %9.4f  %9.4f  %s\n', satname,datestr(t0,0),OE,tt);
   end
   
   % Fill output structure

   ntle=ntle+1;
   tle(ntle).name=deblank(satname);       % Satellite name
   tle(ntle).satid={ satid , ...          % Other satellite identifiers
       classification, deblank(intldesg) }; 
   tle(ntle).ephtype={ ephtype, elnum };  % Ephemeris type and element number
   
   tle(ntle).year=epochyr;                % Epoch Year
   tle(ntle).epoch=epochdays;             % Epoch Day (Day of the year and fractional portion of the day)
   tle(ntle).t0=t0;                       % Epoch (Matlab datenumber)

   tle(ntle).ecc0=ecc0;                   % Eccentricity [-]
   tle(ntle).inc0=inc0*d2r;               % Inclination [rad]
   tle(ntle).raan0=raan0*d2r;             % Right Ascension of the Ascending Node [rad]
   tle(ntle).argp0=argp0*d2r;             % Argument of periapsis [rad]
   
   tle(ntle).m0=m0*d2r;                   % Mean anomaly [rad]
   tle(ntle).n0=n0*2*pi;                  % Mean motion [rad/day]
   tle(ntle).ndot=ndot*2*pi;              % First Time Derivative of the Mean Motion [rad/day^2]
   tle(ntle).nddot=nddot*2*pi;            % Second Time Derivative of Mean Motion  [rad/day^3]
   tle(ntle).bstar=bstar;                 % BSTAR drag term
   tle(ntle).revnum=revnum;               % Revolution number [-]

   tle(ntle).a0=a0*1e3;                   % Semi-major axis [m]    
   tle(ntle).e0=e0*d2r;                   % Eccentric anomaly [rad]
      
end
fclose(fid);

end

function E=eanomaly(M,e)
%EANOMALY    Calculate the eccentric anomaly from mean anomaly
%   E=EANOMALY(M,e)  Calculate the eccentric anomaly E from mean anomaly
%   M and eccentricity e.

err = 1e-10;            %Calculation Error

E0 = M; t =1;
itt = 0;
while(t) 
      E =  M + e*sind(E0);
      if ( abs(E - E0) < err)
          t = 0;
      end
      E0 = E;
      itt = itt+1;
end

end

    
