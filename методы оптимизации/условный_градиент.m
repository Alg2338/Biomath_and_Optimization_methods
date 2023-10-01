%% условный градиент

clc;
f = @(x) (5 * x(2))^2 +(2 * x(1) - 0.7)^2 - 2;
grad = @(x) [8*x(1)-2.8; 50*x(2)-20];
A = [1 0;
     0 1;
     -1 0;
     0 -1];
b = [1; 1; 1; 1];

x = linspace(-10, 10, 100);
y = linspace(-10, 10, 100);

J = @(x, y) (5 * y ).^2 +(2 * x - 0.7).^2 - 2;
[X, Y] = meshgrid(x, y);
contour(X, Y, J(X, Y), 50);
hold on;
%drawSet(@(x) rho(x), 100);
axis equal;
hold on;

u0= [5; 1];
eps = 0.1;
[u_min, J_min] = my_grad(A, b, u0, eps, f, grad)



clear;

