function [tle,fname]=tleget(tleset,tlefile)
%TLEGET   Retrieve NORAD Two Line Elements from www.celestrak.com
%   [TLE,FNAME]=TLEGET(TLESET) download current NORAD two line elements from 
%   Celestrak (www.celestrack.com) and stores the two line elements in
%   the structure array TLE. TLESET is the name of the two-line element 
%   set, or a cell array of strings, with the names of two-line element
%   sets. This can be the name of a file on Celestrak, or, a family of
%   satellites such as GNSS. FNAME is a cell array with the filenames
%   from Celestrak that have been read. 
%
%   The mean orbital elements contained in TLE are Earth-centered-inertial (ECI)
%   coordinates with respect to the true equator of date and the mean equinox 
%   of date. They do not include the effect of nutation. The TLE are only 
%   compatible with the SGP4 and SDP4 orbit propagators of NORAD and should
%   not be used with other orbit propagators.
%
%   [TLE,FNAME]=TLEGET(TLESET,SAVEFILE) also saves the downloaded NORAD two
%   line elements in the file with name SAVENAME. If SAVENAME is not given
%   the two-line elements are written to a temporary file which is read
%   by TLEREAD.
%
%   TLEGET(TLESET,SAVEFILE) saves the downloaded NORAD two line 
%   elements to the file SAVEFILE. SAVEFILE can be read by a seperate
%   call to TLEREAD.
%
%   TLEGET(TLESET) saves the downloaded NORAD two line elements to the 
%   file TLESET.txt in the current directory.
%
%   Examples:
%      tle=tleget('gnss');
%      tleget('gps','gps-20131102.tle')
%
%   Common sets of two line elements are 'GPS' ('gps-ops'), 'GLONASS'
%   ('glo-ops'), 'GALILEO', 'BEIDOU', 'SBAS'; 'GNSS' or 'SATNAV' to do 
%   all satellite navigation systems; 'resource' for Earth resource
%   satellites, etc. For a full list see the Celestrack website.
%   
%   See also TLEREAD, TLEPLOT, TLE2VEC and TLE2AZEL.

%   (c) Hans van der Marel, Delft University of Technology, 2013-2015.
%
%   Created:     2 November 2013 by Hans van der Marel
%   Modified:    3 August 2015 by Hans van der Marel
%                  - improved checking of input / output arguments
%                  - minor correctcions

celestrakurl='http://celestrak.com/NORAD/elements/';

if nargin < 1 || nargin > 2
   disp('TLEGET requires one or two input arguments')
   if nargout > 0, tle=[]; fname={};, end
   return
end

if nargin  == 1
   if nargout == 0
     if ischar(tleset)
       disp('TLEGET with no output arguments and one input argument, save TLE set to default name')
       tlefile=[ tleset '.txt' ];
     else
       disp('TLEGET with no output arguments, and one input argument, input argument must be string')
       return;
     end
     savefile=true;
   else
     savefile=false;
   end
else
   savefile=true;
end

if ischar(tleset)
   tleset=cellstr(tleset);
end
if ~iscellstr(tleset)
   error('Input to TLEGET must be string or a cell array of strings');
end

% Expand family names for the two-line element sets

kk=0;
for k=1:length(tleset)
   tlesetk=lower(tleset{k});
   switch tlesetk
      case {'gnss','satnav'} 
        tleset2{kk+1}='gps-ops';
        tleset2{kk+2}='glo-ops';
        tleset2{kk+3}='galileo';
        tleset2{kk+4}='beidou';
        tleset2{kk+5}='sbas';
        kk=kk+5;
      case {'gps'}
        tleset2{kk+1}='gps-ops';
        kk=kk+1;
      case {'glonass','glo'}
        tleset2{kk+1}='glo-ops';
        kk=kk+1;
      otherwise
        tleset2{kk+1}=tlesetk;
        kk=kk+1;
   end
end

% Get two line elements from Celestrak and store in string s

s=[];

kk=0;
for k=1:length(tleset2)
   tlesetk=tleset2{k};
   fullurl=[celestrakurl tlesetk '.txt'];
   [sk,status]=urlread(fullurl);
   if status
      kk=kk+1;
      fname{kk}=tlesetk;
      s=[s sk];
      disp(['Downloaded ' tlesetk '.txt from ' celestrakurl ]);
   else
      disp(['TLEGET: Warning, could not retrieve ' tlesetk '.txt from ' celestrakurl ]);
   end
end

% Save string to (temporary) file and read using TLEREAD

if ~savefile
   tlefile=tempname;
end
fid=fopen(tlefile,'w');
fprintf(fid,'%s',s);
fclose(fid);

if nargout >= 1
   tle=tleread(tlefile);
end

if savefile
  disp(['Saved TLE to ' tlefile ]); 
else
  delete(tlefile);
end

return
