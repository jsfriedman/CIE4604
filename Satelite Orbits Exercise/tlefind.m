function varargout=tlefind(tle,satid)
%TLEFIND    Find named satellites in the NORAD Two Line Elements.
%   ISAT=TLEFIND(TLE,SATID) finds the satellite with name SATID in structure
%   array TLE with NORAD two line elements, and returns the position(s) in 
%   ISAT. TLE is a structure array that can be read by TLEREAD. SATID must
%   be a character string, cell array of strings, or numeric.
%
%   [ISAT,SATIDS]=TLEFIND(TLE,SATID) also returns the satellite names in the 
%   character cell array SATIDS.
%
%   TLEFIND(TLE,SATID) does the same, but does not return the indices,
%   but instead prints them with the names of the satellites found.
%
%   TLEFIND(TLE) just prints all the satellite names.
%
%   Example:
%
%      tle=tleread('resource.txt');   % read two-line elements
%      tlefind(tle,'RADAR');
%      isat=tlefind(tle,'RADARSAT-2');
%
%   See also TLEREAD, TLEGET, TLE2VEC and TLEPLOT.
%
%   (c) Hans van der Marel, Delft Universtiy of Technology, 2012-2015

%   Created:    30 Dec 2012 by Hans van der Marel
%   Modified:    3 August 2015 by Hans van der Marel
%                  - made this a seperate function (was part of TLE2VEC1)
%                  - enhanced functionality (wildcards/regexp)

% Check the input arguments

if nargin < 1 || nargin > 2
  error('Inproper syntax, use idx=tlefind(tel,satid)')
end
if ~isstruct(tle)
  error('The first argument must be a structure array')
end
if nargin < 2
  satid=[1:length(tle)]';
end

% Find the satellite SATID

if ischar(satid)
  satid=cellstr(satid);
end
if iscellstr(satid)
  satlist=deblank({tle(:).name});
  isat=[];
  for k=1:length(satid)
    satidk=deblank(satid{k});
    % implement a few aliases for otherwise difficult to find families
    switch upper(satidk)
        case 'GLONASS'
           satidk='COSMOS';
        case 'GALILEO'
           satidk='E\d\d';
    end
    if exist('regexpi','builtin')
      satidk=strrep(strrep(satidk,'(','\('),')','\)');
      isati=find(~cellfun('isempty',feval(@regexpi,satlist,satidk,'match','once')));
      isat=[isat isati];
    else
      isat=[isat find(strncmpi(satidk,satlist,length(satidk)))];
    end
  end
elseif isnumeric(satid)
  isat=satid;
else
  error('The second argument must be a character string, cell array of strings, or numeric')
end

if isempty(isat)
  disp(['Satellite ' satid{1} '... not found, please try again.'])
end

% print or return output

if nargout == 0 && length(isat) >= 1
  disp(['Found ' num2str(length(isat)) ' satellites:'])
  for i=1:length(isat)
    disp(['  '  tle(isat(i)).name '   (' num2str(isat(i)) ')'])
  end
  disp('')
end

if nargout == 1
  varargout={isat};
elseif nargout > 1
  satlist=deblank({tle(isat).name})';
  varargout={ isat satlist };
end

end