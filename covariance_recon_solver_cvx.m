function [C,vio] = covariance_recon_solver_cvx(it, Q_list, C_list, v_list, cqi_list, codebook)

% This function is used to estimate the analytic center of problem 1. 

[NA, ~, ~] = size(Q_list); 
[~, M] = size(codebook);

N = norm(Q_list(:,:,1)*codebook(:, v_list(1)));
cvx_begin quiet 
variable C(NA, NA) complex semidefinite
obj = log_det(C) - 0.1*norm_nuc(C);
for i = 1:it
    v_m0 = codebook(:, v_list(i));
    cqi = cqi_list(i);
    Q_t = Q_list(:, :, i);
    if i > 1
        C_t = C_list(:, :, i-1);
    else
        C_t = eye(NA);
    end 
    temp = 0;
    for m = 1:M
        if m ~= v_list(i)
            v_m = codebook(:, m);
           if i == 1 || real(v_m0'*Q_t'*C_t*Q_t*v_m0) - real(v_m'*Q_t'*C_t*Q_t*v_m)<= 1e-10
                  temp=temp+1;

           end
        end 
    end
end 
maximize obj;
subject to

   trace(C) <= 1; 
   lower_bound = -inf;
    for i = 1:it
        v_m0 = codebook(:, v_list(i));
        cqi = cqi_list(i);
        Q_t = Q_list(:, :, i);
        real(v_m0'*Q_t'*C*Q_t*v_m0) == cqi;
        if cqi/real(trace(Q_t*(v_m0*v_m0')*Q_t'))>lower_bound
            lower_bound = cqi/real(trace(Q_t*(v_m0*v_m0')*Q_t'));
        end
    end 
   trace(C) >= lower_bound;

cvx_end


vio = 0;
for i =1:it
         Q_t = Q_list(:, :, i);
         v_m0 = codebook(:, v_list(i));
        for m = 1:M
            v_m = codebook(:, m);
            if (v_m'*Q_t'*C*Q_t*v_m - v_m0'*Q_t'*C*Q_t*v_m0)>eps
                vio = vio+1;
            end
        end
end

end 