%% Demonstration script to simulate a random walk process
%
% Compute random walk variable x from recursive formula
%
%    x(k) = x(k-1)+w(k)
%
% with w(k) normal distributed with expectation 0 and variance qk.
% The variance qk is an input parameter.
% 
% The length of the series is n, and the number of simulations is m.
%
% We assume that the time interval dt is one, then qk=q*dt, with q the 
% power spectral dimension.

qk=1;      % This should be one of the input parameters
n=2000;    % Should also be input
m=100;

%% Simulate matrix with values (three different methods)
%
% Each realization is in a column of x. The size of x is [n,m].

tic;
x=nan(n,m);      % Preallocate memory 
w=sqrt(qk)*randn(n,m);  %  N(0,qk) qk is variance;

% for p=1:m
%    x(1,p)=0;   % Initialize
%    for k=2:n
%      x(k,p)=x(k-1,p)+w(k,p);
%    end
% end

x(1,:)=zeros(1,m);
for k=2:n
   x(k,:)=x(k-1,:)+w(k,:);
end

% x=cumsum(w);

toc

% Plot the realizations

figure
plot(x)

%% Compute the emperical mean and standard deviation

xmean=mean(x,2);
xstd=std(x,0,2);

hold on
h(1)=plot(xmean,'b','linewidth',2);
h(2)=plot(xmean+2*xstd,'r','linewidth',2);
plot(xmean-2*xstd,'r','linewidth',2);

%% compute the formal variance and add to plot

sx=zeros(n,1);
for k=2:n
  sx(k)=sx(k-1)+qk;
end
sx=sqrt(sx);

h(3)=plot(2*sx,'k','linewidth',2);
plot(-2*sx,'k','linewidth',2);

xlabel('time [s]')
ylabel('rw [-]')
title([ 'Random-walk process (#sim=',num2str(m),', q_k=',num2str(qk),')'])

legend(h,'Mean','Emp.std.','Std')