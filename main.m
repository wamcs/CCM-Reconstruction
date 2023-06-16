function [ratioc,rmses,vios,C_list] = main(C,T)
% @ parameter
% C: the ground-truth channel covariance matrix
% T: the number of communication rounds.


%---------Parameter Setting---------
NA = 32; NP = 8; max_it =  T; NU = 2; 


%---------Codebook---------
codebook = get_codebook(NP); 

%---------C---------

C = C ./ trace(C); 
[u1, d1]=eigs(C,1);


%%---------Q---------
Q_list = zeros(NA, NP, max_it);
Q_list(:, :, 1) = init_Q(C); 

v_list = zeros(max_it, 1); % feedback list 
cqi_list = zeros(max_it, 1); 

rmses = [];  % rmses lists 
ratioc = [];


C_list = zeros(NA, NA, max_it); 
vios = [];


for it = 1 : max_it  
    Q_t = Q_list(:,:,it);
    [v_list(it), cqi_list(it)] = get_feedback(Q_t, C, codebook); % feedback
    fprintf('the feedback is %d \n', v_list(it));
    
    [C_est,vio] = covariance_recon_solver_cvx(it, Q_list, C_list, v_list, cqi_list, codebook);
    vios = [vios,vio];
    C_list(:, :, it) = C_est;  
    
    [w,d_est] = eigs(C_est,1);
    ratio =  real(w'*C*w/(d1));
    ratioc = [ratioc, ratio];
    rmse = norm(C_est - C, 'fro')/(NA^2);
    rmses = [rmses, rmse]; 
    if rmse <= 1e-4
        break; 
    end 
    fprintf('the ratio is %f at iteration %d \n', ratio, it);
    fprintf('the rmse is %f at iteration %d \n', rmse, it);
    
    Q_list  (:,:,it+1) = pilot_design(it, C_est, Q_list, codebook); % problem p2

end
end
