function [u_min, J_min] = simplex(A, b, c, v)
    J_min = -inf;
    u_min=[];
    r=rank(A);
    eps = 0.001;
    i = 0;
    while true
        i = i+1;
        [Jb, ~] = find(v >= eps);
        [Jf, ~] = find(v < eps);
        if numel(Jb) < r
            m = r - numel(Jb);
            Jb = [Jb; Jf(1:m)];
            Jf = Jf(m+1:end);
        end
        B = A(1:end, Jb);
        F = A(1:end, Jf);
        vb = b;
        cb = c(Jb);
        cf = c(Jf);
        vb = v(Jb);
        C = B \ F;
        delta = C' * cb - cf;
        Jf_plus = find(delta > eps); 
        if isempty(Jf_plus)
            break;
        end
        k = min(Jf_plus);
        gamma = C(1:end,k);
        Ik = find(gamma > eps);
        if isempty(Ik)
            u_min = [];
            J_min = -inf;
            return;
        end
        theta = min(vb(Ik) ./ gamma(Ik));
        u_min = zeros(size(v));
        u_min(Jb) = vb - theta * gamma;
        u_min(Jf(k)) = theta;
        J_min = u_min' * c;
        v = u_min;               
    end
    disp(['Algorithm converged in ', num2str(i), ' iteration(s)']);
end
