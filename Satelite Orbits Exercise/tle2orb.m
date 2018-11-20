function orb=tle2orb(tle,t,propagation)
%TLE2ORB    Satellite orbital elements from NORAD Two Line Elements.
%   ORB=TLE2ORB(TLE,T,PROPAGATION) computes the satellite orbital elements
%   ORB at time T from NORAD two line elements in the TLE structure.
%   The ouput ORB is a matrix with on each row the six orbital elements
%   corresponding to time T: semi-major axis [m], eccentricity [-], 
%   inclination [rad], right ascension of ascending node [rad], argument of 
%   periapsis [rad] and true anonaly [rad]. TLE is an element of the 
%   structure array that can be read by TLEREAD.  T must be a vector with 
%   Matlab date numbers in UT1. PROPAGATION is the orbit propagation method
%   that is used to compute ORB. ORB can be used by ORB2VEC from the
%   CRSUTIL toolbox to compute the state vector (position and velocity).
%
%   Valid PROPAGATION methods are:
%
%      J2     Include secular effects of J2, but nothing else. (default)
%             This gives acceptable results for plotting purposes and 
%             planning computations. 
%
%      NOJ2   Ignores the effect of J2 on the orbit propagation.
%             This method should only be used for educational purposes.
%
%      SGP4   SGP4 orbit propagator of NORAD (not yet supported). 
%             Formally, TLE are only compatible with the SGP4 orbit 
%             propagator of NORAD, but this is a complicated computation
%             that has not yet been implemented.
%
%   ORB=TLE2ORB(TLE,T) uses the default propagation method J2.
%
%   See also TLEGET, TLEREAD, TLEFIND, TLEPLOT, KEPLERM, TLE2VEC and ORB2VEC.
%
%   (c) Delft University of Technology, 2017

%   Created:    13 September 2017 by Hans van der Marel
%   Modified:   13 September 2017 by Hans van der Marel
%                  - split of from TLE2VEC
%                  - include J2 in orbit propagation

% Constants

J2 = 0.00108262998905;     % J2 
Re = 6378136;              % [m]   radius of the Earth
mu = 3986004418e5;         % [m^3/s^2] gravitational constant of the Earth

% Check the input arguments

if nargin < 2 || nargin > 3
  error('Inproper syntax, use orb=tle2vec(tle,t[,propagation])')
end
if ~isstruct(tle)
  error('The first argument must be a structure array')
end
if nargin < 3
  propagation='J2';
end

t=tledatenum(t);

% Initialize orbit propagation

nepoch=length(t);
doy=floor(tle.epoch);
t0=datenum(tle.year,1,doy)+tle.epoch-doy;
    
% Orbit propagation

switch upper(propagation)
   case 'J2'
       % Compute rate of change of orbital elements
       %
       %   draan/dt = s.*cos(inclination)
       %   dargp/dt = -0.5*s.*(5*cos(inclination-1).^2)
       %   dM/dt = -0.5*s.*sqrt(1-e.^2).*(3*cos(inclination).^2 -1)
       %
       % with s=-J2*3/2*sqrt(mu/a^3)*(Re/p)^2
       %
       % dM/dt is not needed for two line element propagation, but computed nevertheless. 

       p=tle.a0*(1-tle.ecc0^2);
       s=-J2*3/2*sqrt(mu/tle.a0.^3)*(Re/p)^2;  % identical to s=-J2*3/2*n0*(Re/p)^2;

       odot = s * cos(tle.inc0) * 86400;
       wdot = -0.5 * s * ( 5*cos(tle.inc0)^2 -1 ) * 86400;
       mdot = -0.5 * s * sqrt(1-tle.ecc0^2) * ( 3*cos(tle.inc0)^2 -1 ) * 86400;

       raan=tle.raan0+odot*(t-t0);                 % Right ascension of ascending node [rad]
       argp=tle.argp0+wdot*(t-t0);                 % Argument of periapsis [rad]

       m=tle.m0+tle.n0*(t-t0);                     % Mean anomaly [rad]
       [~,nu]=keplerm(m+argp,tle.ecc0);            % Solve Kepler's equation for E + argp
       nu=nu-argp;                                 % True anomaly  [rad]

       orb= [ ...                                  % Matrix of orbit elements
           repmat([ tle.a0 tle.ecc0 tle.inc0 ],[nepoch 1]) ...
           raan argp nu ];            
       
   case 'NOJ2'
       % Very simple orbit propagation ignoring effect of J2 (use with
       % extreme caution, only for educational purposes)
       m=tle.m0+tle.n0*(t-t0);                     % Mean anomaly 
       [~,nu]=keplerm(m,tle.ecc0);                 % True anomaly
       orb= [ ...                                  % Matrix of orbit elements
           repmat([ tle.a0 tle.ecc0 tle.inc0 tle.raan0 tle.argp0 ],[nepoch 1]) ...
           nu ];            

   case 'SGP4'
       error('SGP4 propagation model is not yet supported')
       
   otherwise
       error('Illegal propagation method')

end


return

