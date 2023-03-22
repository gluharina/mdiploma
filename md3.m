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

for A=0.1:0.1:100
    count = count + 1;
    [t,yp] = ode45(@(t,y) func2(t,y,p,T,A), [0 N], [1; 1; 0]);
%     psi = zeros(1,length(M)); u = zeros(1,length(M));
%     ksi = sin(t(:,1));
%     uu = yp(:,1) - (a*yp(:,1) - yp(:,1).*yp(:,2)/(1+aa*yp(:,1)) - e*yp(:,1).*yp(:,1)) - ksi;
    sum = 0;
    for i = 1500:1700
        s(i-1499) = yp(i,1) - p;
    end
    sum = var(s);
    sumc(count) = sum;
    ac(count) = A;
%     fprintf('%.5f',sum);
%    plot(t, yp(:,1) - p,'Linewidth',3);
%     hold on;
end

plot(ac, sumc,'Linewidth',3);
xlabel("амплитуда гармонического шума"),ylabel("дисперсия по макропеременной");
% % 
% figure;
% plot(t, yp(:,1)/n,'Linewidth',3);
% title('\psi(t) ode45');

function out = func2(t,y,p,T,A)
    a = 0.8; aa = 0.8; e = 0.8; g = 0.2; b = 0.05;
    n = 1; k = 1;

%     A = 1;
    w = 10;
    ksi = A*sin(t*w);

    psi = y(1) - p;
    psiend = psi + k.*y(3);
    f1 = a*y(1) - y(1).*y(2)/(1+aa*y(1)) - e*y(1).*y(1);
    u = y(3).*(k*k*n -1) - psiend./T - f1;
    y1 = f1 + u + ksi;
    y2 = -g*y(2) + y(1).*y(2)/(1+aa*y(1)) - b*y(2).*y(2);
    y3 = n*psi;
       
    out = [y1; y2; y3;];
end
