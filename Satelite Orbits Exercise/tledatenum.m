function t=tledatenum(range)
%TLEDATENUM   Compute Matlab datenumbers from a date range.
%   T=TLEDATENUM(DATERANGE) computes the Matlab datenumbers from the
%   three element cell array DATERANGE with  start date, end date (or 
%   duration) and data interval. The duration and data interval are in 
%   minutes.  DATERANGE can also be a character or cell array with dates.
%
%   Example:
%      t=tledatenum({'2013-9-13 0:00', 24*60 ,1}) 
%
%   See also TLE2VEC, TLE2VEC1, TLEPLOT and TLEPLOT1.
%
%   (c) Delft University of Technology, 2012-2015

%   Created:    30 Dec 2012 by Hans van der Marel
%   Modified:    4 August 2015 by Hans van der Marel
%                  - split of seperate function TLEDATENUM from TLEPLOT1

if iscell(range) && length(range) == 3 && ~iscellstr(range)

  daterange(1)=datenum(range{1}); 
  if ischar(range{2})
     daterange(2)=datenum(range{2}); 
  else
     daterange(2)=daterange(1)+range{2}/(24*60);
  end
  daterange(3)=range{3}/(24*60);

  fprintf('range:    %s - %s\n',datestr(daterange(1),0),datestr(daterange(2),0));
  fprintf('stepsize: %s\n',datestr(daterange(3),13));

  t=(daterange(1):daterange(3):daterange(2))';

elseif ischar(range) || iscellstr(range)
  
  t=datenum(range);

elseif isnumeric(range)

  t=range;

else
  
  error('Invalid input argument')
  
end

if size(t,2) > 1
  t=t';
end
if size(t,2) > 1
  error('The input argument must be a scalar or vector')
end

end