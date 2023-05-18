clc, clear, 
close all

T = 1;
a = .32; aa = (1-0.21)/1; g = .60; b = 0.0011; e = 0.0087;
n = 1; k = 1; 

y = [1689,24,0];
N = 100;
h = 0.1;
M = 0:h:N;

const = 0; A = 0; w = 0; sigma = 0; mu = 0;

for goal = 1:2
    if goal == 1 
        strg = 'x_1* = '; 
    else
        strg = 'p = ';
    end
    
    for pr = 1:14
   
        if pr == 1
            const = 10;
            mode = 1;
            str = 'const = 10';
        elseif pr == 2
            const = 1000000;
            mode = 1;
            str = 'const = 1000000 ';
        elseif pr == 3
            const = 0.000001;
            mode = 1;
            str = 'const = 0.000001';
        elseif pr == 4
            mode = 2;
            A = 0.00001;
            w = 0.1;
            str = 'A = 0.00001 w = 0.1 ';
        elseif pr == 5
            mode = 2;
            A = 10;
            w = 0.00001;
            str = 'A = 10 w = 0.00001';
        elseif pr == 6
            mode = 2;
            A = 10;
            w = 0.1;
            str = 'A = 10 w = 0.1';
        elseif pr == 7
            mode = 2;
            A = 1000;
            w = 100000;
            str = 'A = 1000 w = 100000';
        elseif pr == 8
            mode = 2;
            A = 10000000;
            w = 1;
            str = 'A = 10000000 w = 1';
        elseif pr == 9
            mode = 3;
            mu = 0;
            sigma = 0.1;
            str = '\mu = 0 \sigma = 0.1';
        elseif pr == 10
            mode = 3;
            mu = 0;
            sigma =1; str = '\mu = 0 \sigma = 1';
        elseif pr == 11
            mode = 3;
            mu = 0;
            sigma =10; str = '\mu = 0 \sigma = 10';
        elseif pr == 12
            mode = 3;
            mu = 0;
            sigma =100; str = '\mu = 0 \sigma = 100 ';
        elseif pr == 13
            mode = 3;
            mu = 0;
            sigma =1000;   str = '\mu = 0 \sigma = 1000 ';
        elseif pr == 14
            mode = 2;
            A = 10000;
            w = 10;
            str = 'A = 10000 w = 10';
        end

        [t,yp] = ode45(@(t,y) func2(t,y,T,mode,goal,const,A,w,mu,sigma), [0 N], [1689; 24; 0]);

        if mode == 1
            ksi = const;
        elseif mode == 2
            ksi = A*sin(t(:,1)*w);
        else
            ksi = normrnd(mu,sigma);
        end
        
         if goal == 1
           p = 10000; 
           psi = yp(:,1) - p;
        else
            p = 10;
            psi = yp(:,1) - p*yp(:,2);
         end
        
        uu = yp(:,1) - (a*yp(:,1) - yp(:,1).*yp(:,2)./(1+aa*yp(:,1)) - e*yp(:,1).*yp(:,1)) - ksi;

        figure;
        subplot(2,2,2);
        plot(t, uu(:,1),'Linewidth',3);
        title('Управление'); legend('u(t)');
        ylabel('Численность популяций');
        xlabel('Время, t');
        % saveas(gcf,'C:\Users\categ\Desktop\диплом\коэф2\цель1\ksi = 10000\управление.png');

        subplot(2,2,4);
        plot(t, psi/n,'Linewidth',3);title('\psi(t)'); legend('\psi(t)');
        xlabel('Время, t');title('\psi(t)'); ylabel('Численность популяций');
        % saveas(gcf,'C:\Users\categ\Desktop\диплом\коэф2\цель1\ksi = 10000\цель.png');

        subplot(2,2,[3,1]);
        plot(t,yp,'Linewidth',3);
        xlabel('Время, t'); ylabel('Численность популяций');
        lgd = legend({'x_{1}', 'x_{2}','z'});
        title(lgd,[str ' ' strg  num2str(p)])
        title('Графики популяций');
        % saveas(gcf,'C:\Users\categ\Desktop\диплом\коэф2\цель1\ksi = 10000\график.png');

        if goal == 1 && mode == 1
            sgtitle('Графики по цели x_1 -> x_1* с внешним шумом - постоянной');
        elseif goal == 1 && mode == 2
            sgtitle('Графики по цели x_1 -> x_1* с внешним шумом A(sin(t)*w)');
        elseif goal == 1 && mode == 3
            sgtitle('Графики по цели x_1 -> x_1* с нормально распределенным внешним шумом');
        elseif goal == 2 && mode == 1
            sgtitle('Графики по цели x_1 -> px_2 с внешним шумом - постоянной');
        elseif goal == 2 && mode == 2
            sgtitle('Графики по цели x_1 -> px_2 с внешним шумом A(sin(t)*w)');
        elseif goal == 2 && mode == 3
            sgtitle('Графики по цели x_1 -> px_2 с нормально распределенным внешним шумом');
        end
    end
end


function out = func2(t,y,T,mode,goal,const,A,w,mu,sigma)
    a = .32; aa = (1-0.21)/1;  g =.60; b = 0.0011; e = 0.0087;
    n = 1; k = 1;
    
    if mode == 1
        ksi = const;
    elseif mode == 2 
        ksi = A*sin(t*w);
    elseif mode == 3
        ksi = normrnd(mu,sigma);
    end
    
    if goal == 1
       p = 10000; 
       psi = y(1) - p;
    elseif goal == 2
        p = 10;
        psi = y(1) - p*y(2);
    end
    
    psiend = psi + k.*y(3);
    f1 = a*y(1) - y(1).*y(2)/(1+aa*y(1)) - e*y(1).*y(1);
    u = y(3).*(k*k*n -1) - psiend./T - f1;
    y1 = f1 + u + ksi;
    y2 = -g*y(2) + y(1).*y(2)/(1+aa*y(1)) - b*y(2).*y(2);
    y3 = n*psi;
    out = [y1; y2; y3;];
end
