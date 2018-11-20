function [xsat,vsat]=tle2vec1(tle,t,satid,propagation)
%TLE2VEC1    Satellite position and velocity from NORAD Two Line Elements.
%   [XSAT,VSAT]=TLE2VEC1(TLE,T,SATID) computes the satellite position XSAT
%   and velocity VSAT from NORAD two line elements for satellite SATID at times
%   T. TLE is a structure array that can be read by TLEREAD. T must be a vector 
%   with Matlab date numbers in UT1, a character string or cell array of 
%   strings.
%
%   The dimension of XSAT and VSAT is 2. The first dimenison is time T,
%   the second dimension is the coordinates. Therefore, size(XSAT) and size(VSAT)
%   is [ m 3], with m=lenght(T).
%
%   [...]=TLE2VEC1(....,PROPAGATION) selects an alternative propagation 
%   method PROPAGATION. The default propagation method is J2 which
%   includes the secular effects of J2. This gives acceptable results
%   for plotting purposes and planning computations. Formally, TLE are only 
%   compatible with the SGP4 orbit propagator of NORAD, but this is not
%   (yet) supported by this function. The other available option is NOJ2, 
%   which ignores the effect of J2 on the orbit propagation and should only
%   be used for educational purposes.
%
%   Example:
%
%      tle=tleread('resource.txt');   % read two-line elements
%      tle2vec1(tle,[],'R');          % satellite names starting with R
%      [xsat,vsat]=tle2vec1(tle,'2013-9-13 0:00','RADARSAT-2') 
%
%      t0=datenum('2013-9-13 0:00');  % vector with date numbers ...
%      t=[t0:1/(24*60):t0+1];
%      [xsat,vsat]=tle2vec1(tle,t,'RADARSAT-2') 
%
%   See also TLEREAD, TLEGET, TLE2VEC, KEPLERM, TLE2ORB and ORB2VEC.
%
%   (c) Hans van der Marel, Delft Universtiy of Technology, 2012-2017

%   Created:    30 Dec 2012 by Hans van der Marel
%   Modified:   13 September 2013 by Hans van der Marel
%                  - select unique SATID
%                  - added remark about J2 effects
%                3 August 2015 by Hans van der Marel
%                  - renamed to tle2vec1
%                  - added example
%                  - minor improvements, more flexible inputs
%                  - call to TLEFIND and TLEDATENUM
%               13 Sep 2017 by Hans van der Marel
%                  - change of tle structure format
%                  - call to TLE2ORB for orbit propagation
%                  - include propagation method as fourth argument

% Check the input arguments

if nargin < 3 || nargin > 4
  error('Inproper syntax, use [xsat,vsat]=tle2vec1(tel,t,satid)')
end

if ~isstruct(tle)
  error('The first argument must be a structure array')
end

if nargin < 4
  propagation='J2';
end

t=tledatenum(t);

% Find the satellite SATID

isat=tlefind(tle,satid);
if isempty(isat)
  %disp(['Satellite ' satid ' not found, please try again.'])
  if nargout > 0, xsat=[];vsat=[];,end
  return
elseif length(isat) > 1
  tlefind(tle,satid);
  disp(['Satellite ' satid ' not unique, please try again.'])
  if nargout > 0, xsat=[];vsat=[];,end
  return
end

% Compute satellite state vectors (position and velocity in ECI)

orb=tle2orb(tle(isat),t,propagation);      % Matrix of orbit elements
vec=orb2vec(orb);                          % State vector

xsat=vec(:,1:3);
vsat=vec(:,4:6);

return