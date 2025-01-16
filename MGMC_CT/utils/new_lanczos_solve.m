function [x, iter] = new_lanczos_solve(A_phi1,B, b, epsilon, iter_max,level)
%
    iter = 0;
    n = size(b, 1);
    x = zeros(n,1);
    w = b;
    r_norm = norm(w);
    v = w / r_norm;
    beta = 0;
    p = zeros(n, 1);
%     T2=A_phi1';
T2=A_phi1{2};
T_A_phi1=A_phi1{1};
% T=T_A_phi1*T2;
    while r_norm >= epsilon && iter < iter_max
        iter = iter + 1;
        if iter == 1
            T1=T2*v;
            TAV=T_A_phi1*T1;           
            w =TAV+B.*v;
%              w=T_A_phi1*(T2*v)+B.*v;
            lambda = 0;
            xi = r_norm;
        else
            T1=T2*v;
            TAV=T_A_phi1*T1;
%             w =T_A_phi1*(T2*v)+B.*v - beta * v_prev;
              w =TAV+B.*v - beta * v_prev;
            lambda = beta / u;            % u == u_{k-1,k-1}
            xi = - lambda * xi;
        end
        v_prev = v;
        alpha = dot(w, v);                % Modified Gram-Schmit
        u = alpha - lambda * beta;        % u := u_{k,k}
        p = (v - beta * p) / u;
        x = x + xi * p;
        w = w - alpha * v;
        beta = norm(w);                   % beta := beta_{k+1}
        v = w / beta;                     % v := v_{k+1}
        r_norm = abs(beta * xi / u);
    end
end