clc, clear, 
close all

p = 10;
T = 3;
a = 0.8; aa = 0.8; e = 0.8; g = 0.2; b = 0.05;
n = 1; k = 1; 

y = [1,1,0];
N = 100;
h = 0.01;
M = 0:h:N;
count = 0;

[t1,yp1] = ode45(@(t,y) func1(y,p,T,1), [0 N], [1; 1; 0]);
s2 = std(yp1(:,1));

% for sigma=0.001:0.001:10
    count = count + 1;
    [t,yp] = ode45(@(t,y) func2(y,p,T,sigma), [0 N], [1; 1; 0]);
%     plot(t, yp);
%     hold on;
%     psi = zeros(1,length(M)); u = zeros(1,length(M));
%     ksi = sin(t(:,1));
%     uu = yp(:,1) - (a*yp(:,1) - yp(:,1).*yp(:,2)/(1+aa*yp(:,1)) - e*yp(:,1).*yp(:,1)) - ksi;

    s1 = std(yp(:,1));
    s11(count) = s1;
    sums = s1/s2; 
    sigmap(count) = sums;
    
    sum = 0;
    for i = 1500:1700
        s(i-1499) = yp(i,1)-p;
    end
    sum = var(s);
    sumc(count) = sum;
   
    fprintf('%i ',  count)
%     
%     if count > 100
%         break
%     end
end
plot(sigmap, sumc,'.');
% hist(sigmap, sumc);
xlabel("отношение сигнал-шум"),ylabel("разброс по макропеременной");

function out = func2(y,p,T,sigma)
    a = 0.8; aa = 0.8; e = 0.8; g = 0.2; b = 0.05;
    n = 1; k = 1;

    mu = 0;
%     sigma =100;
    ksi = normrnd(mu,sigma);
    
    psi = y(1) - p;
    psiend = psi + k.*y(3);
    f1 = a*y(1) - y(1).*y(2)/(1+aa*y(1)) - e*y(1).*y(1);
    u = y(3).*(k*k*n -1) - psiend./T - f1;
    y1 = f1 + u + ksi;
    y2 = -g*y(2) + y(1).*y(2)/(1+aa*y(1)) - b*y(2).*y(2);
    y3 = n*psi;
       
    out = [y1; y2; y3;];
end

function out = func1(y,p,T,sigma)
    a = 0.8; aa = 0.8; e = 0.8; g = 0.2; b = 0.05;
    n = 1; k = 1;

    mu = 0;
%     sigma =100;
    ksi = normrnd(mu,sigma);
    ksi=0;
    
    psi = y(1) - p;
    psiend = psi + k.*y(3);
    f1 = a*y(1) - y(1).*y(2)/(1+aa*y(1)) - e*y(1).*y(1);
    u = y(3).*(k*k*n -1) - psiend./T - f1;
    y1 = f1 + u + ksi;
    y2 = -g*y(2) + y(1).*y(2)/(1+aa*y(1)) - b*y(2).*y(2);
    y3 = n*psi;
       
    out = [y1; y2; y3;];
end
