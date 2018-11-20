% Jonathan Friedman 4799380
% CIE4604 Assignment 1
% test
dt = 1;
q = 0.5;
n = 2000; % length of walk
m = 100; % number of walks
x = zeros(m,n); % walk results

qk = q*dt;
for j = 1:m
    for i = 2:n
        x(j,i) = x(j,i-1) + sqrt(qk)*randn();
    end
end