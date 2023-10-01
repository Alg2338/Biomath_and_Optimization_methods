function [u_min, J_min] = my_grad(A, b, u0, eps, J, grad)
    us = [u0];
    alphas = linspace(0, 1, 100);
    nmax = 1000;
    u_prev = u0;
    der = grad(u_prev);
    lin = @(x) der(1) * x(1) + der(2) * x(2);
    y = fmincon(lin, u_prev, A, b);
    g = @(a) J(u_prev + a * (y - u_prev));
    Js = arrayfun(g, alphas);
    [~, al_opt] = min(Js);
    u_curr = u_prev + alphas(al_opt) * (y - u_prev)
    us = [us, u_curr];
    n = 1;
    while (n < nmax) && (norm(u_curr - u_prev) > eps)
        u_prev = u_curr;
        der = grad(u_prev);
        lin = @(x) der(1) * x(1) + der(2) * x(2);
        y = fmincon(lin, u_prev, A, b);
        g = @(a) J(u_prev + a * (y - u_prev));
        Js = arrayfun(g, alphas);
        [~, al_opt] = min(Js);
        u_curr = u_prev + alphas(al_opt) * (y - u_prev)
        us = [us, u_curr];
        n = n + 1;
    end
    
    us
    plot(us(1,:), us(2, :), '-o')
    u_min = u_curr;
    J_min = J(u_curr);
end

