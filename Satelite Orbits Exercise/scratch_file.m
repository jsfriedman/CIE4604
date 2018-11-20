qk = 0.5;
n = 2e3;
m = 100;
x = nan(n,m);
w = sqrt(qk)*randn(n,m);
        
% for p=1:m
%     x(1,p) = 0;
%     for k=2:n
%         x(k,p) = x(k-1,p)+w(k,p);
%     end
% end

x(1,:) = zeros(1,m);
for k=2:n
    x(k,:) = x(k-1,:)+w(k,:);
end
    
plot(x);
hold on

% compute formal variance
vx = zeros(n,1);
for k=2:n
    vx(k) = vx(k-1)+qk;
end

xmean = mean(x,2);
xstd = std(x,0,2);
plot(xmean,'k', 'linewidth',3);
plot(1.96*xstd+xmean,'r','linewidth',3);
plot(-1.96*xstd-xmean, 'r', 'linewidth',3);
plot(1.96*sqrt(vx), 'k', 'linewidth',2.5);



qk = 0.5;
n = 2000;
m = 100;
x = cumsum(randn(m,n)*sqrt(qk));
std = sqrt(cumsum(qk*ones(1,n)));

plot(x);
hold on
