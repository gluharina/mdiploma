clear all;

syms T;

str = input('введите символы: ','s');
C = strsplit(str);
i=0;
while i<size(C,2)
    i = i+1;
    syms(C(i));
end

n = input ('колво уравнений в системе ');

z = zeros(n,3);
z = sym(z);
i=0;
while i<n
    i = i+1;
    z(i,1) = str2sym(strcat('x',string(i),'(t)'));
    z(i,2) = input (strcat('x',string(i),' = '));
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

syms T;
if k > 0 
    sol = subs(diff(psi), old, new);
    EL = simplify(T*sol + psi);
    pretty(solve(EL,u))
else 
    syms fi T1 T2;
    i=0;
    psi_1 = psi-psi;
    while i<n
        i = i+1;
        if strfind(string(psi), string(z(i,1))) > 0 
%            psi_1 = psi_1 - subs(fi, t, z(i,1));
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
    disp('-------------------->');
    disp('u = ');
    pretty(solve(EL_1,u))
    sol_2 = subs(diff(psi), old, new);
    sol_3 = subs(sol_2, x_u, fi);
    EL_2 = (T2*sol_3 + psi);
    fii = solve(EL_2,fi,'ReturnConditions',true);
    disp('-------------------->');
    disp('fi = '); pretty(fii.fi)
    
    diff_fi = subs(fi*diff(x_fi), old, new);
    disp('-------------------->');
    disp('diff_fi = ');
    pretty(diff_fi)
end


