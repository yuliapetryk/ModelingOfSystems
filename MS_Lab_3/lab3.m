clear
clc

#Часовий інтервал
dt = 0.2;
t0 = 0;
tk = 50;
eps = 1e-6;

#Вказуємо відомі параметри або наближення 
c = [0.1, 0.3, 0.1, 0.12]';
m = [12, 19, 18]';

# Матриця утворена функцією чутливості
function result = matrixA(m, c)
     A = [
        0, 1, 0, 0, 0, 0;
        -(c(2) + c(1)) / m(1), 0, c(2) / m(1), 0, 0, 0;
        0, 0, 0, 1, 0, 0;
        c(2) / m(2), 0, -(c(2) + c(3)) / m(2), 0, c(3) / m(2), 0;
        0, 0, 0, 0, 0, 1;
        0, 0, c(3) / m(3), 0, -(c(4) + c(3)) / m(3), 0
    ];
  result = A;
end

#Метод Рунге-Кутти для функції чутливості
function result = findU (A, U, y, m, c, dt)
  D = calculateModelDerivatives(y, m, c);
  k1 = dt * (A * U + D);
  k2 = dt * (A * (U + k1/2) + D);
  k3 = dt * (A * (U + k2/2) + D);
  k4 = dt * (A * (U + k3) + D);

  result = U + (k1 + 2*k2 + 2*k3 + k4)/6;
end

#Обчислюємо похідні для моделі
function result = calculateModelDerivatives(y, m, c)
  db = zeros(6, 3);
  db(2, 1) = -(y(1))/m(1);
  db(4, 2) = (y(5)-y(3))/m(2);
  db(6, 2) = (y(3)-y(5))/m(3);
  db(4, 3) =-(c(2)*y(1)-c(2)*y(3)-c(3)*y(3)+c(3)*y(5))/m(2)^2;

  result = db;
end

#Зчитуємо дані з файлу 
data = dlmread('y8.txt');
[size_m, size_n] = size(data);

#Задаємо показник якості ідентивікації 
I = inf;

i = 1;

while I > eps
  y = data(:, 1);
  dy = zeros(6, 1);
  U = zeros(6, 3);
  integral1 = zeros(3, 3);
  integral2 = zeros(3, 1);
  I = 0.0;
  
  A = matrixA(m, c);

  for j = 2:size_n
      
    #Метод Рунге-Кутти для знаходження y(t, beta0)
    k1 = dt*(A * y);
    k2 = dt*(A * (y + k1/2));
    k3 = dt*(A *(y + k2/2));
    k4 = dt*(A * (y + k3));
    
    y_new = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    
    U_new = findU(A, U, y_new, m, c, dt);
 
    dy_new = data(:, j) - y_new;
       
    integral1 = integral1 + dt*(U'*U + U_new'*U_new) / 2;
    integral2 = integral2 + dt*(U'*dy + U_new'*dy_new) / 2;

    #Обраховуємо значення показника якості ідентивікації 
    I = I + dt*(dy'*dy + dy_new'*dy_new) / 2;

    U = U_new; 
    y = y_new; 
    dy = dy_new;
    
  end    
  
  dBeta = pinv(integral1) * integral2;
  c(1) = c(1) + dBeta(1);
  c(3) = c(3) + dBeta(2);
  m(2) = m(2) + dBeta(3);
  
  i = i + 1;
end

#Виводимо значення показника якості ідентивікації на поточному кроці 
  disp('Показник якостi iдентифiкацiї параметрiв:')
  disp(['I = ', num2str(I)])
  
#Виводимо шукані параметри
disp('Шукані параметри:');
for i = 1:4
    fprintf('  c%d = %s\n', i, num2str(c(i)));
end
disp('');
for i = 1:3
    fprintf('  m%d = %s\n', i, num2str(m(i)));
end

