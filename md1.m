clc, clear, close all

p = 10;
T = 3;
a = 0.8; aa = 0.8; e = 0.8; g = 0.2; b = 0.05;
n = 1; k = 1; 

y = [1,1,0];
N = 100;
h = 0.01;
M = 0:h:N;
[t,yp] = ode45(@(t,y) func2(t,y,p,T), [0 N], [1; 1; 0]);
psi = zeros(1,length(M)); u = zeros(1,length(M));

for i=1:length(M)-1
    ksi = 2;
 %ksi = 100000000*sin(i*h*0.1 + 3);
%     ksi = randn*10;
    psi(i) = y(i,1) - p;
    psiend = psi(i) + k.*y(i,3);
    f1 = a.*y(i,1) - y(i,1).*y(i,2)/(1+aa.*y(i,1)) - e.*y(i,1).*y(i,1);
    u(i) = y(i,3).*(k*k*n -1) - psiend./T - f1;
    y1 = f1 + u(i) + ksi;
    y2 = -g*y(i,2) + y(i,1).*y(i,2)/(1+aa.*y(i,1)) - b*y(i,2)*y(i,2);
    y3 = n*psi(i);
    
    y(i+1,1) = y(i,1) + h*y1; 
    y(i+1,2) = y(i,2) + h*y2;
    y(i+1,3) = y(i,3) + h*y3;
end

plotting2(M,y,psi,u,t,yp);

function out = func2(t,y,p,T)
    a = 0.8; aa = 0.8; e = 0.8; g = 0.2; b = 0.05;
    n = 1; k = 1; 
%     ksi = 2;
    ksi = sin(t)/5;
%     ksi = randn;
    psi = y(1) - p;
    psiend = psi + k.*y(3);
    f1 = a*y(1) - y(1).*y(2)/(1+aa*y(1)) - e*y(1).*y(1);
    u = y(3).*(k*k*n -1) - psiend./T - f1;
    y1 = f1 + u + ksi;
    y2 = -g*y(2) + y(1).*y(2)/(1+aa*y(1)) - b*y(2).*y(2);
    y3 = n*psi;
    out = [y1; y2; y3;];
end

function plotting2(M,y,psi,u,t,yp)
    figure;
    plot(t,yp,'Linewidth',3);
    xlabel('t');
    legend({'Y_{1}', 'Y_{2}','Z'});
    title('1');

    figure;
    plot(M(1:end-1),psi(1:end-1),'Linewidth',3);
    title('\psi(t)');

    figure;
    plot(M(1:end-1),u(1:end-1),'Linewidth',3);
    title('Управление');
    
    figure;
    plot(M(1:end-1),y(1:end-1,1),M(1:end-1),y(1:end-1,3),'Linewidth',3);
    title('\xi and z');
    legend({'f_{1} + u + (\xi)','z'});
    
    figure;
    plot(M,y,'Linewidth',3);
    xlabel('t');
    legend({'y_{1}', 'y_{2}','z'});
end
