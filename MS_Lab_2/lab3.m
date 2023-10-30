clear, clc, close all 
 
% Matrix norm
function result = norm(matrix)
    result = sum(sum(abs(matrix)));
end
 
% Moore-Penrose
function result = MoorePenrose(matrix, eps)
    
    delta = 1;
    eps_now = eps+1;
    E = eye(size(matrix, 1));
    pseudo_inv= transpose(matrix) * inv(matrix * transpose(matrix) + delta * E);
    while eps_now > eps
        new_pseudo_inv= pseudo_inv;
        delta = delta/2;
        pseudo_inv = transpose(matrix) * inv(matrix * transpose(matrix) + delta*delta * E);
        eps_now = norm(pseudo_inv - new_pseudo_inv);
    end
  
    result = pseudo_inv;
end
 
% Greville
function result = Greville(matrix)
   
    pseudo_inv = 0;
    
    a = matrix(1,:)';
   
    matrix_now = a';
  
    if(a'*a == 0)
       pseudo_inv  = a;
    else
        pseudo_inv  = a / a' * a;
    end
   
    i = 2;
    while i <= size(matrix, 1)
       
        Z = eye(size(pseudo_inv , 1))-pseudo_inv * matrix_now;
        R = pseudo_inv *pseudo_inv ';
        row = matrix(i,:)';
        matrix_now = [ matrix_now ; row'];
        aZa = row'*Z*row;
        aRa = 1+row'*R*row;
        if(aZa == 0)
            pseudo_inv  = [(pseudo_inv -(R*row*row'*pseudo_inv )/aRa),(R*row)/aRa];
        else
            pseudo_inv = [(pseudo_inv -(Z*row*row'*pseudo_inv )/aZa),(Z*row)/aZa];
        end
        i = i + 1;
    end
    result =pseudo_inv ;
end

%Pseudoinversion check
function checkProperties(pseudoInverse, originalMatrix)
    %AA+A = A;
    firstProperty = isequal(round(originalMatrix * pseudoInverse * originalMatrix), 
                            round(originalMatrix));
                            
    %A+AA+ = A+;
    secondProperty = isequal(round(pseudoInverse * originalMatrix * pseudoInverse),
                             round(pseudoInverse));
                             
    %AA+ – симетрична матриця розмiрностi m × m
    thirdProperty = isequal(round(originalMatrix * pseudoInverse), 
                            round(transpose(originalMatrix * pseudoInverse)));
                            
    %A+A – симетрична матриця розмiрностi n × n
    fourthProperty = isequal(round(pseudoInverse * originalMatrix),
                             round(transpose(pseudoInverse * originalMatrix)));
                             
    fprintf('Properties check results %d %d %d %d\n', firstProperty, secondProperty, thirdProperty,  fourthProperty);  
end

%Initial signals
X = double(imread('x1.bmp'));
Y = double(imread('y8.bmp'));
 
%Pseudo-inverse matrices
resultMoore = MoorePenrose(X, 0.00001);
resultGreville = Greville(X);

%Check properties
checkProperties(resultMoore, X)
checkProperties(resultGreville, X)

%Random matrix
V = rand(size(Y,1), size(X,1));

%Сalculation of the linear operator
A1 = Y*resultMoore + V * (eye(size(X, 1))-X*resultMoore)';
A2 = Y*resultGreville + V * (eye(size(X, 1))-X*resultGreville)';
 
%Display the results
figure;
plot(1,1);
imshow(uint8(X));
title('X1');

figure;
plot(1,4);
imshow(uint8(Y));
title('Y8');

figure;
plot(1,4);
imshow(uint8(A1*X));
title('Moore-Penrose');

figure;
plot(1,4);
imshow(uint8(A2*X));
title('Greville');
