%% Newton Test 1
f2 = @(x) (x.^4) / 4 + 2 * x.^2 - 3 * x;
grad2 = @(x) x.^3 + 4 * x - 3;
hess2 = @(x) 3 * x.^2 + 4;
x0 = 1;
x_min = newton_min(f2, grad2, hess2, x0);
x_min = x_min(1, :);
disp('Newtons solution:');
disp(x_min); 
t = linspace(0, 1.5, 30);
plot(t, f2(t), x_min, f2(x_min), 'o');

%% Newton Test 2
f3 = @(x) (x(1).^4) / 4 - (x(1).^2) / 2 + x(2).^2;
grad3 = @(x) [x(1).^3 - x(1); 2 * x(2)];
hess3 = @(x) [3 * x(1).^2 - 1, 0; 0, 2];
x0 = [3; 4];
x_min = newton_min(f3, grad3, hess3, x0);
sz = 300;
other_f3 = @(x, y) x.^4 / 4 - x.^2 / 2 + y .^ 2;
[X, Y] = meshgrid(linspace(-6, 6, sz));
plot(x_min(1, :), x_min(2, :), 'o-r');
hold on;
contour(X, Y, other_f3(X, Y), other_f3(x_min(1,:), x_min(2, :)), '--');
hold on;
axis([-4 4 -6 6]);
axis equal;
disp('Newtons solution:');
disp(x_min);


%% simplex
clc;
clear;
A = [1 0 0 1 -1;1 1 0 2 0;0 0 1 1 0];
b = [1; 3; 1];
c = [1;1;-1;0;1];
z0 = [3; 0; 1;0;2];
[u, J] = simplex(A, b, c, z0);
disp("J = ");
disp(J);
disp("u = ");
disp(u);


