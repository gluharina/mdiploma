clear all;
syms a aa e y p b T xx x1(t) x2(t) u;
syms vm ks m d g rm ly T xx u x1(t) x2(t) x3(t) p
% f1 = -vm*x1/(ks + x1)*x2 + m*x2 + d*x3 + (1 - g)*x3*rm.*(1-exp(T)^(-ly*x2)) + u; %нутриенты
% f2 = vm*x1/(ks + x1)*x2 -  m*x2 - x3*rm*(1-exp(T)^(-ly*x2)) ; %фитоплантон 
% f3 = rm*g*x3*(1-exp(T)^(-ly*x2)) - d*x3; %зоопланктон 
% f1 = a*x1 - (x1*x2)/(1 + aa*x1) - e*x1^2 + u;
% f2 = -y*x2 + x1*x2/(1 + aa*x1) - b*x2^2;

n = input ('колво уравнений в системе ');

z = zeros(n,3);
z = sym(z);
i=0;
while i<n
    i = i+1;
    z(i,1) = str2sym(strcat('x',string(i),'(t)'));
    z(i,2) = input (strcat('x',string(i),' =  '));
    if strfind(string(z(i,2)), 'u') > 0
        z(i,3) = 1;
    end
end


psi = input ('цель = ');
psi_z = zeros(n,1);
i=0;k = 0;
while i<n
    i = i+1;
    if strfind(string(psi), string(z(i,1))) > 0 
        if z(i,3) > 0
            k = k + 1;
        end
    end
end 

i=0;
old = sym(zeros(1,2));
new = sym(zeros(1,2));
while i<n
    i = i+1;
    old(i) = diff(z(i,1));
    new(i) = z(i,2);
end
    
if k > 0 
    sol = subs(diff(psi), old, new);
    EL = simplify(T*sol + psi);
    pretty(solve(EL,u))
else 
    syms fi(t) T1 T2;
    i=0;
    psi_1 = psi-psi;
    while i<n
        i = i+1;
        if strfind(string(psi), string(z(i,1))) > 0 
            psi_1 = psi_1 - fi;
            x_fi = z(i,1);
        end
        if z(i,3) > 0
           psi_1 = psi_1 + z(i,1);
           x_u = z(i,1);
        end
    end    
    sol_1 = subs(diff(psi_1), old, new);
    EL_1 = simplify(T1*sol_1 + psi_1);
    pretty(solve(EL_1,u))
    sol_2 = subs(diff(psi), old, new);
    sol_3 = subs(sol_2, x_u, fi);
    EL_2 = simplify(T2*sol_3 + psi);
    
    solve(EL_2,fi)
    
    diff_fi = subs(fi*diff(x_fi), old, new);
    pretty(diff_fi)
end
