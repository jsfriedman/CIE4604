function tleplot1(tle,range,satid,objcrd)
%TLEPLOT1   Plot satellite position and velocity from NORAD Two Line Elements.
%   TLEPLOT1(TLE,RANGE,SATID,OBJCRD) plots the satellite position and velocities
%   in various ways. TLE is a structure array with two line elements read by 
%   TLEREAD. RANGE is a cell array with the start date, end date (or duration) 
%   and data interval. The duration and data interval are in minutes. SATID is 
%   a character string with the name of the satellite to plot. OBJCRD are the 
%   geographical coordinates with latitude and longitude (in degrees), and 
%   height (in meters) of the observer or object on Earth.
%
%   Example:
%      tle=tleread('resource.txt')
%      tleplot1(tle,{'2013-9-13 0:00', 24*60 ,1},'RADARSAT-2',[ 52 4.8  0 ]) 
%
%   Files with TLE's can be obtained from www.celestrak.com
%
%   See also TLEGET, TLEREAD, TLEFIND, TLEDATENUM, TLE2VEC1, KEPLERM and ORB2VEC.
%
%   (c) Delft University of Technology, 2012-2015

%   Created:    30 Dec 2012 by Hans van der Marel
%   Modified:   13 September 2013 by Hans van der Marel
%                  - adapted from rxtle2vec for single satellites for CIE4606
%                    course and named TLEPLOT
%                4 August 2015 by Hans van der Marel
%                  - renamed to TLEPLOT1 for single satellites, TLEPLOT is
%                    doing multiple satellites
%                  - call to TLE2VEC1 instead of multi satellite TLE2VEC 

% Constants

Re = 6378136;         % [m]   radius of the Earth   
mu = 3986004418e5;    % [m^3/s^2] gravitational constant of the Earth
Me = 7.2921151467e-5; % [rad/sec] rotational velocity of the Earth
c  = 299792458;       % [m/s] speed of light        

% Compute the epoch times t from the cell array with start, end and data interval

t=tledatenum(range);
nepoch=length(t);

% ----------------------------------------------------------------------------
% Compute and plot satellite position and velocity in ECI
% ----------------------------------------------------------------------------

% Compute satellite state vectors (position and velocity in ECI) using tle2vec.m

[xsat,vsat]=tle2vec1(tle,t,satid);
if isempty(xsat), return;, end

% Compute the satellite radius and velocity in an ECI

rsat=sqrt(sum(xsat.^2,2));
velsat=sqrt(sum(vsat.^2,2));

% Plot satellite position and velocity in ECI (in km and km/s))

figure('Name','ECI Position','NumberTitle','off')
plot(t,xsat./1000)
hold on
plot(t,rsat./1000,'k','LineWidth',2)
title([ satid ' position in ECI'])
ylabel('[km]')
datetick('x')
legend('X','Y','Z','r')

figure('Name','ECI Velocity','NumberTitle','off')
plot(t,vsat./1000)
title([ satid ' velocity in ECI'])
hold on
plot(t,velsat./1000,'k','LineWidth',2)
ylabel('[km/s]')
datetick('x')
legend('Vx','Vy','Vz','v')

% ----------------------------------------------------------------------------
% Compute position and velocity of the observer in ECI
% ----------------------------------------------------------------------------

% The position of the observer (latitude, longitude and height) is given in ECEF
% using the input array objcrd

lat=objcrd(1)*pi/180;     % convert latitude from degrees to radians
lon=objcrd(2)*pi/180;     % convert longitude from degrees to radians
Rs=Re+objcrd(3);          % convert height to radius (from CoM) using mean Earth radius Re 

% Position of the observer in ECEF (assume latitude and longitude are for spherical Earth)

xobjECEF=Rs* [ cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
vobjECEF=[ 0 0 0 ];
%
% The transformation from an ECEF to ECI is a simple rotation around the z-axis
% a. the rotation angle is GMST (Greenwhich Mean Stellar Time) minus UT1
% b. the rotation around the z-axis can be implemented by replacing the
%    longitude (in ECEF) by local stellar time (lst) in the ECI
%
% The times given in t are in UTC, which is close to UT1 (max 0.9 s difference),
% which is not important for a plotting application

% Compute GMST from UT1, for the first epoch in t, using the Matlab function 
% ut2gmst. The second output returned by ut2gmst is the rotational velocity
% omegae of the Earth in rev/day

[gst0,omegae]=ut2gmst(t(1));

% Compute local stellar time (in radians) from the longitude, GMST at the initial
% epoch and the rotational velocity of the Earth (times elapsed time). Note that
% lst is an array, while lon is a scalar)

lst=lon+gst0+2*pi*omegae*(t-t(1));

% Compute position and velocity of the observer in ECI using lst (position and 
% velocity in an ECI change all the time, unlike in a ECEF)

xobj=NaN(nepoch,3);      % pre-allocate memory, makes Matlab faster
vobj=NaN(nepoch,3);

xobj(:,1)=Rs*cos(lat)*cos(lst);
xobj(:,2)=Rs*cos(lat)*sin(lst);
xobj(:,3)=Rs*sin(lat)*ones(nepoch,1);
vobj(:,1)=-Rs*cos(lat)*Me*sin(lst);
vobj(:,2)=Rs*cos(lat)*Me*cos(lst);
vobj(:,3)=zeros(nepoch,1);

% Plot observer position (Not really interesting, but easy to do, use previous
% plots as example)


% ----------------------------------------------------------------------------
% Satellite position and velocity from the observer (object) point of view
% ----------------------------------------------------------------------------

% position and velocity vectors, range and range rates, from object to satellites

xobj2sat=xsat-xobj;
vobj2sat=vsat-vobj;

robj2sat=sqrt(sum(xobj2sat.^2,2));
rrobj2sat=dot(vobj2sat,xobj2sat./repmat(robj2sat,[1 3]),2);

% Note that the range rate rrobj2sat is not the same as the relative velocity
% sqrt(sum(vobj2sat.^2,2)), these are different things

% normal vector (vertical) and unit direction vector to satellite from observer

robj=sqrt(sum(xobj.^2,2));               % range to the object (observer)
n0=xobj./repmat(robj,[1 3]);             % normal vector from object (observer)
ers=xobj2sat./repmat(robj2sat,[1 3]);    % init direction vector from observer to satellite

% zenith angle and azimuth of satellite (as seen from object wrt to radial direction)

ip=dot(n0,ers,2);
zenith=acos(ip);
azi= atan2( -n0(:,2).*ers(:,1) + n0(:,1).*ers(:,2) ,   ip .* -n0(:,3)  + ers(:,3) ) ;
azi=mod(azi+2*pi,2*pi);

% Elevation angle and satellite visibility (if elevation angle > 0))

elevation=pi/2-zenith;

cutoff=0;
visible=elevation > cutoff;

% Plot elevation, azimuth, range and range rate

figure('Name','Viewing angles...','NumberTitle','off')
subplot(2,2,1)
plot(t,elevation.*180/pi,'b:')
hold on
plot(t(visible),elevation(visible).*180/pi,'b.')
ylabel('Elevation angle [deg]')
datetick('x')
title(satid)

subplot(2,2,3)
plot(t,azi.*180/pi,'b:')
hold on
plot(t(visible),azi(visible).*180/pi,'b.')
ylabel('Azimuth angle [deg]')
datetick('x')
title(satid)

subplot(2,2,2)
plot(t,robj2sat./1000,'b:')
hold on
plot(t(visible),robj2sat(visible)./1000,'b.')
ylabel('Range [km]')
datetick('x')
title(satid)

subplot(2,2,4)
plot(t,rrobj2sat./1000,'b:')
hold on
plot(t(visible),rrobj2sat(visible)./1000,'b.')
ylabel('Range rate [km/s]')
datetick('x')
title(satid)

% Make a skyplot (polarplot of azimuth and zenith angle)

skyplot(t,azi,zenith,cutoff*180/pi);     % Internal function

% ----------------------------------------------------------------------------
% Plot satellite orbit tracks in ECI (right-ascension and declination)
% ----------------------------------------------------------------------------

% Compute right-ascension and declination of the observer and satellite

aobj=atan2(xobj(:,2),xobj(:,1))*180/pi;
dobj=atand(xobj(:,3)./sqrt(xobj(:,2).^2+xobj(:,1).^2));
asat=atan2(xsat(:,2),xsat(:,1))*180/pi;
dsat=atand(xsat(:,3)./sqrt(xsat(:,2).^2+xsat(:,1).^2));

% Plot right ascension and declination 

figure('Name','ECI tracks','NumberTitle','off');
plot(aobj,dobj,'g.','MarkerSize',8)
hold on
plot(asat,dsat,'.','MarkerSize',8)
plot(aobj(visible),dobj(visible),'g+','MarkerSize',8)
plot(asat(visible),dsat(visible),'+','MarkerSize',8)
axis([ -180 +180 -90 90]);
legend('observer',satid)
xlabel('Right ascension')
ylabel('Declination')
title([ satid ' orbit track in ECI'])


% ----------------------------------------------------------------------------
% Satellite ground tracks
% ----------------------------------------------------------------------------

% Substract GMST from right-ascension of observer and object in ECI to get the 
% longitude in ECEF

[gst0,omegae]=ut2gmst(t(1));
gst=gst0*180/pi+360*omegae*(t-t(1));

lobj=aobj-gst-360*round((aobj-gst)./360);   % must be in the range [-180,+180]
lsat=asat-gst-360*round((asat-gst)./360);

figure('Name','Ground Tracks','NumberTitle','off')
plot(lobj,dobj,'g*','MarkerSize',8)
hold on;
axis([ -180 +180 -90 90]);
if exist('coast.mat')==2
   coast=load('coast');
   line(coast.long,coast.lat,'color',[.35,.35,.35]);
end
plot(lsat,dsat,'.','MarkerSize',8);
plot(lsat(visible),dsat(visible),'+','MarkerSize',8);
legend('observer',satid)
xlabel('Longitude')
ylabel('Latitude')
title([ satid ' ground tracks'])

% ----------------------------------------------------------------------------
% 3D plot (ECI)
% ----------------------------------------------------------------------------

f = figure('Name','3D-orbit (ECI)','NumberTitle','off');
plot_orbit_3D(xsat,xobj)
legend('observer',satid)
title('3D satellite orbit (ECI)')

% ----------------------------------------------------------------------------
% 3D plot (ECEF)
% ----------------------------------------------------------------------------

% Compute rotation angle (GMST) around Z-axis

[gst0,omegae]=ut2gmst(t(1));
gst=gst0+2*pi*omegae*(t-t(1));
lst=lon+gst;

% Rotate satellite positions round z-axis (ECI->ECEF))

xsate(:,1)=  cos(gst).*xsat(:,1) + sin(gst).*xsat(:,2);
xsate(:,2)= -sin(gst).*xsat(:,1) + cos(gst).*xsat(:,2);
xsate(:,3)=  xsat(:,3);

% Rotate observer positions round z-axis (ECI->ECEF) 

xobje(:,1)=  cos(gst).*xobj(:,1) + sin(gst).*xobj(:,2);
xobje(:,2)= -sin(gst).*xobj(:,1) + cos(gst).*xobj(:,2);
xobje(:,3)=  xobj(:,3);

f = figure('Name','3D-orbit (ECEF)','NumberTitle','off');
plot_orbit_3D(xsate,xobje)
legend('observer',satid)
title('3D satellite orbit (ECEF)')
plot3(xobje(:,1)/1000,xobje(:,2)/1000,xobje(:,3)/1000,'r*')  % to improve visibility in plot


end

% -----------------------------------------------------------------------------
% Internal functions for the more complicated plots
% -----------------------------------------------------------------------------

function h=skyplot(epoch,azi,zen,cutoff)
%SKYPLOT   Create skyplot
%   SKYPLOT(epoch,azi,zen,cutoff) creates a polar plot with the elevation and 
%   azimuth of a satellite. AZI is an array with the azimuth (radians), ZEN 
%   and array with the zenith angle (radians), and CUTOFF is the cutoff 
%   elevation (degrees)


xx = [-1 0 -1].';
yy = [.4 0 -.4].';
arrow = xx + yy.*sqrt(-1);

h1 = figure('Units','normalized', ...
    'Position',[0.2 0.3 0.6 0.6], ...
    'NumberTitle','off',...
    'Name','Skyplot', ...
    'Color',[1 1 1]);

h=gca;
set(h,'Units','normalized', ...
	  'Position',[0.05 0.05 0.90 0.90], ...
	  'Color',[1 1 1], ...
	  'XColor',[1 1 1], ...
	  'YColor',[1 1 1], ...
	  'Xgrid','off', ...
	  'Ygrid','off', ...
	  'Box','off');

%  az = [0:1:360];
%  el = 90*ones(size(az));
%  [x,y] = pol2cart (deg2rad(az),el);
%  patch('xdata',x,'ydata',y, ...
%             'edgecolor',[0 0 0],'facecolor',[1 1 1],...
%             'handlevisibility','off','parent',h);

axis square
hold on;

% -----------------------------------------
% --- Limits, grid, labels and such ... ---
% -----------------------------------------

set (gca,'Xlim',[-90 90]);
set (gca,'Ylim',[-90 90]);

set (gca,'FontWeight','bold');
set (gca,'FontSize',16);
set (gca,'LineWidth',2);

for i = 0:30:330;
  [x,y] = pol2cart (-deg2rad(i-90),90);
  plot ([0 x], [0 y],'k-');
  [x,y] = pol2cart (-deg2rad(i-92),94);
  h = text (x,y,num2str(i));
  set (h,'Color',[0 0 0]);    
  set (h,'FontSize',12);
  set (h,'FontWeight','bold');
  set (h,'Rotation',-i);
end;

for i = [0:15:90 cutoff];
  az = [0:1:360];
  el = (90-i)*ones(size(az));
  [x,y] = pol2cart (deg2rad(az),el);
  h = plot (x,y,'k-');
  set (h,'LineWidth',1);
  if i == cutoff;
    set (h,'LineStyle',':');
  else
    if i ~= 0;
      h = text (2,90-i+3,num2str(i));
      set (h,'Color',[0 0 0]);    
      set (h,'FontSize',12);
      set (h,'FontWeight','bold');
    end;
  end;
end;

% -------------------
% --- Actual plot ---
% -------------------

lcol = get(gca,'ColorOrder');

dt=min(diff(epoch));

idx1 = find( zen < (90-cutoff)*pi/180);
  
if ~isempty(idx1);
    
    idx2 = [0 ;find(diff(epoch(idx1)) > dt*3 ); length(idx1) ];
    
    for j = 1:length(idx2)-1

      idx3=idx1(idx2(j)+1:idx2(j+1)) ;
            
      [x,y] = pol2cart(-azi(idx3)+pi/2, 180/pi*zen(idx3));
      h = plot (x,y);
      set (h,'LineStyle','-');
      set (h,'LineWidth',3);
      set (h,'Color',lcol(mod(i,size(lcol,1))+1,:));
      
      [tx,ty] = pol2cart(-azi(idx3(end))+pi/2, 180/pi*zen(idx3(end)));
      if length(idx3) > 1
        [txx,tyy] = pol2cart(-azi(idx3(end-1))+pi/2, 180/pi*zen(idx3(end-1)));
        dd=sqrt((tx-txx)^2+(ty-tyy)^2);
        z = ((tx-txx)/dd + ((ty-tyy)/dd).*sqrt(-1)).';
        a = arrow * z;
        h=plot(tx+3*real(a),ty+3*imag(a));
        set (h,'LineStyle','-');
        set (h,'LineWidth',3);
        set (h,'Color',lcol(mod(i,size(lcol,1))+1,:));
        tx=tx+6*(tx-txx)/dd;
        ty=ty+5*(ty-tyy)/dd;
      else
        tx = tx + 5; ty = ty + 2;
      end
%      h  = text(tx,ty,num2str(prns(i)));
%      set (h,'VerticalAlignment','Middle','HorizontalAlignment','Center')
%      set (h,'Color',lcol(mod(i,size(lcol,1))+1,:));
%      set (h,'FontSize',18);
%      set (h,'FontWeight','bold');
    end;
end

end

function f=plot_orbit_3D(xsat,xobj)
%PLOT_ORBIT_3D   Plots satellite orbits in 3D around the Earth.

Re = 6378136;         % [m]   radius of the Earth   

[X,Y,Z] = sphere(20);
X=Re*X/1000;
Y=Re*Y/1000;
Z=Re*Z/1000;

colormap('Lines')
hold on
plot3(xobj(:,1)/1000,xobj(:,2)/1000,xobj(:,3)/1000,'r','LineWidth',2)
plot3(xsat(:,1)/1000,xsat(:,2)/1000,xsat(:,3)/1000) % 'b')
mesh(X,Y,Z,2*ones(size(X)),'CDataMapping','direct')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
axis equal
view([45 45]);

end
