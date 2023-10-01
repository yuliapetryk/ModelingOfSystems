clear %

#Load data from file
y = load('f8.txt');

#Initialize value
dt = 0.01; 
T = 5;
t=0:dt:T;
df=1/T;
N = length(y);
n = length(t);
f = 0:df:round(n/2) * df;

#Graph of the initial signal
figure
plot(t,y), grid

#Graph of the Fourier transform
figure
yf = fft(y)/N;
plot(abs(yf)), grid

#Local extrema
figure 
f = 0:df:round(n/2) * df;
plot(f,abs(yf(1:round(n/2)+1)))

#Fourier function
yf = zeros(1, N);
for m = 1:N
  for j = 1:N
    yf(m) = yf(m) + 1/N*y(j)*exp(1)^(-1i*2*pi/N*m*j);
  end
end
yf = abs(yf);

#Calculate extremums
counter = 0;
extrema = zeros(2,1);
for j = 3:round(N/2)-1
  if (yf(j) > yf(j+1) && yf(j) > yf(j-1) && abs(yf(j)-yf(j+1)) > 1)
    counter = counter + 1;
    extrema(counter) = j*df;
  end
end

#Solve the system of equations
function_sin = sin(2*pi*extrema(1)*t);
a = [sum(t.^6), sum(t.^5), sum(t.^4), sum(function_sin.*t.^3), sum(t.^3);
     sum(t.^5), sum(t.^4), sum(t.^3), sum(function_sin.*t.^2), sum(t.^2);
     sum(t.^4), sum(t.^3), sum(t.^2), sum(function_sin.*t),    sum(t);
     sum(function_sin.*t.^3), sum(function_sin.*t.^2), sum(function_sin.*t),  sum(function_sin.*function_sin), sum(N*function_sin);
     sum(t.^3), sum(t.^2), sum(t), sum(N*function_sin), N];

c = [sum(y.*t.^3), sum(y.*t.^2), sum(y.*t), sum(y.*function_sin),  sum(y)];
a = inv(a)*c'
temp = a'

#Graph of the approximating function
function_aprox = a(1).*t.^3 + a(2).*t.^2 + a(3).*t + a(4).*function_sin +a(5);
figure
plot(t, function_aprox), grid

#Calculate error
error_value = sum((function_aprox-y).^2)


