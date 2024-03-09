clearvars; %очищение переменных

%ввод символов, которые будут использоваться для символического задания системы уравнений
er = 0;
while er < 1
    str = input('введите символы: ','s'); %ввод символов одной строкой через пробел
    C = strsplit(str); %разделение строки на массив подстрок 
    er = input ('подтвердить ввод (1 - да, 0 - нет): '); %проверить ранне введеные символы 
end 

%перевод подстрок в символы
i=0;
while i<size(C,2)
    i = i+1;
    syms(C(i)); %каждый элемент массива подстрок C переводится в sym 
end

n = input ('введите количество уравнений в системе: '); % n - количество уравнений в системе
z = sym(zeros(n,3));% z = sym(z); %массив, в котором содержатся: 1 столбец - левая сторона уравнений, 2 столбец - правая сторона уравнений, 3 столбец - наличие в уравнении управления u
syms T u; %?????????? символы T и u
er = 0;
while er < 1
    disp('введите уравнения (управление u должно быть указано)');
    i=0; u_counter = 0;
    while i <= n
        i = i+1; 
        if i < n + 1
            z(i,1) = str2sym(strcat('x',string(i),'(t)')); %заполняется первый столбец массива z - левая сторона уравнения
            z(i,2) = input (strcat('x',string(i),' = ')); %заполняется второй столбец массива z - правая сторона уравнения
            if strfind(string(z(i,2)), 'u') > 0
                z(i,3) = 1; %заполняется третий столбец массива z - наличие в уравнении управления u
                u_counter = u_counter + 1;
            end
        end
        if u_counter < 1 && i == n %проверка на наличие u в системе уравнений
            disp('управление не было введено');
            i = 0;
        end
    end
    er = input ('подтвердить ввод (1 - да, 0 - нет): ');
end

er = 0;
while er < 1
    psi = input ('цель = '); % psi - цель
    er = input ('подтвердить ввод (1 - да, 0 - нет): ');
end

i = 0; k = 0; %счетчик k = 1 - цель задана по x_i по которому введено управление, k = 0 - цель и управление заданы по разным x_i
while i < n
    i = i+1;
    if strfind(string(psi), string(z(i,1))) > 0 % если в цели присутвует x_i  
        if z(i,3) > 0 %если в управление в системе по x_i
            k = k + 1; 
        end
    end
end 

%вспомогательные массивы для замен
old = sym(zeros(1,2));
new = sym(zeros(1,2));
i = 0;
while i<n
    i = i+1;
    old(i) = diff(z(i,1)); %массив заполняется симовольными значениями типа diff(x_i(t), t)
    new(i) = z(i,2); %массив заполняется правой частью x_i(t)
end

if k > 0 %если счетчик k = 1, то цель задана по x_i по которому введено управление
    sol = subs(diff(psi), old, new); %замена в diff(psi) - diff(x_i(t), t) на правую часть x_i(t)
    EL = simplify(T*sol + psi); %уравнение эйдера-лагранжа
    disp('-------------------->'); %вывод управления
    disp('u = ');
    pretty(solve(EL,u))
else %если счетчик k = 0 - цель и управление заданы по разным x_i
    syms fi T1 T2;
    i=0;
    psi_1 = sym('0');
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
