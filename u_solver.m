clearvars; %очищение переменных

%ввод символов, которые будут использоваться для символического задания системы уравнений
er = 0; %счетчик для повторного ввода 
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
syms T u; 
er = 0; %счетчик для повторного ввода 
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

er = 0;%счетчик для повторного ввода 
while er < 1
    psi = input ('цель = '); % psi - цель
    er = input ('подтвердить ввод (1 - да, 0 - нет): ');
end

i = 0; k = 0; %счетчик k > 0 - цель задана по x_i по которому введено управление, k = 0 - цель и управление заданы по разным x_i
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
    i=0;
    psi_1 = sym('0'); %задаем новую цель_1 = (x_i по которому управление) - fi
    while i < n
        i = i+1;
        if z(i,3) > 0
            psi_1 = psi_1 + z(i,1); %задаем (x_i по которому управление) для новой цели_1
            x_u = z(i,1); %определяем x_i по которому управление
        end
    end
    
    i = 0; m = 0;  
    while i <= n
        i = i+1; 
        if i < n + 1 
            if strfind(string(psi), strcat('x',string(i))) > 0 % проверка на нахождение x_i в первоначальной цели
                if strfind(string(z(i,2)),erase(string(x_u),'(t)')) > 0  % проверка на присутствие в правой части x_i в переменной по которой задано управление
                    m = m + 1; %счетчик m, если m = 0, то понадобится дополнительная макропеременная для передачи управления
                else 
                    x_fi_start = z(i,1); %x_i, по которой задана первоначальная цель 
                end
            end
            if strfind(string(z(i,2)),erase(string(x_u),'(t)')) > 0 % проверка на присутствие в правой части x_i в переменной по которой задано управление
                if string(z(i,1)) ~= string(x_u) % x_i не является переменной по которой задано управление
                    x_fi = z(i,1); %x_i на которую передается управление, если будет использоваться дополнительная макропеременная
                end    
            end    
        end
    end
    if m > 0 % если m > 0 то управление можно вывести с помощью одной дополнительной макроперменной
        syms fi(t) T1 T2 fi_output;
        
        psi_1 = psi_1 - fi; %задаем fi для новой цели_1  
        sol_1 = subs(diff(psi_1), old, new); %замена в diff(psi) - diff(x_i(t), t) на правую часть x_i(t)
        EL_1 = simplify(T1*sol_1 + psi_1); %первое уравнение эйдера-лагранжа для поиска управдения
        disp('-------------------->');
        disp('u = ');
        pretty(solve(EL_1,u)) %выражаем управдение содержащее fi и diff(fi)
        
        sol_2 = subs(diff(psi), old, new);
        sol_3 = subs(sol_2, x_u, fi); %замена (x_i по которому управление) на fi
        EL_2 = (T2*sol_3 + psi); %второе уравнение эйдера-лагранжа для поиска fi 
        EL_2 = subs(EL_2,fi,fi_output); % замена функции fi(t) на переменную fi_output для вывода ее значения
        fii = solve(EL_2,fi_output,'ReturnConditions',true); %выражаем fi
        disp('-------------------->');
        disp('fi = '); pretty(fii.fi_output) %выводим fi без условий

        diff_fi = subs(fi*diff(x_fi), old, new); %выражаем diff(fi)
        disp('-------------------->');
        disp('diff(fi) = ');
        pretty(diff_fi) %выводим diff(fi)
    else % если m = 0 то для вывод управления понадобится больше одной передачи управления черерз макропеременные 
        syms fi_1(t) fi_2(t) T1 T2 T3 fi1_out fi2_out; 
             
        psi_1 = psi_1 - fi_1; %задаем fi_1 для дополнительной цели_1  
        sol_1 = subs(diff(psi_1), old, new); %замена в diff(psi) - diff(x_i(t), t) на правую часть x_i(t)
        EL_1 = simplify(T1*sol_1 + psi_1); %первое уравнение эйдера-лагранжа для поиска управления
        disp('-------------------->');
        disp('u = ');
        pretty(solve(EL_1,u)) %выражаем управдение содержащее fi и diff(fi)
        
        psi_2 = x_fi - fi_2; %задаем fi_2 для дополнительной цели_2  
        sol_2 = subs(diff(psi_2), old, new); %замена в diff(psi) - diff(x_i(t), t) на правую часть x_i(t)
        sol_3 = subs(sol_2, x_u, fi_1);
        EL_2 = (T2*sol_3 + psi_2); %второе уравнение эйдера-лагранжа для поиска управления
        EL_2 = subs(EL_2,fi_1,fi1_out); % замена функции fi_1(t) на переменную fi1_out для вывода ее значения
        fii = solve(EL_2,fi1_out,'ReturnConditions',true); %выражаем fi_1
        disp('-------------------->');
        disp('fi_1 = '); pretty(fii.fi1_out) %выводим fi_1 без условий
                
        sol_4 = subs(diff(psi), old, new); %замена в diff(psi) - diff(x_i(t), t) на правую часть x_i(t)
        sol_5 = subs(sol_4, x_fi, fi_2);
        EL_3 = (T3*sol_5 + psi); %третье уравнение эйдера-лагранжа для поиска управления
        EL_3 = subs(EL_3,fi_2,fi2_out); % замена функции fi_2(t) на переменную fi2_out для вывода ее значения
        fii1 = solve(EL_3,fi2_out,'ReturnConditions',true); %выражаем fi_2
        disp('-------------------->');
        disp('fi_2 = '); pretty(fii1.fi2_out) %выводим fi_2 без условий
        
        diff_fi_1 = subs(fi_1*diff(x_fi), old, new); %выражаем diff(fi_1)
        disp('-------------------->');
        disp('diff(fi_1) = ');
        pretty(diff_fi_1) %выводим diff(fi_1)
        
        diff_fi_2 = subs(fi_2*diff(x_fi_start), old, new); %выражаем diff(fi_2)
        disp('-------------------->');
        disp('diff(fi_2) = ');
        pretty(diff_fi_2) %выводим diff(fi_2)
    end
end
