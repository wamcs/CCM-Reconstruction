function [type_i,type_ii] = benchmark(H,stream)
% C: ground truth 
 
C = H'*H;
 
[~, d]=eigs(C,stream);
d = trace(d);

%% Type I ratio
Q0 = init_Q(C); 
codebook = get_codebook(8); 
[v_idx, ~] = get_feedback(Q0, C, codebook); 
v = codebook(:,v_idx);
% for i = 1:stream
%     v(i) = v(i)/norm(v(i));
% end
% 
% 
% Q0 = Q0/norm(Q0);
%  
baseline_init = real(trace(v'*Q0'*C*Q0*v)/d/(norm(Q0*v))^2);
type_i = baseline_init; 
    
    %% Type II ratio
WW = TypeII_rank1(H, H);
WW = WW/norm(WW); 
type_ii = real(trace(WW'*C*WW)/d); 

% for n = 1:500
%     H = ds(:, :, n, 1) + 1i*ds(:, :, n, 2);
%     C = H'*H;
%     [~, d1]=eigs(C,1);
%     
%     %% Type I ratio
%     Q0 = init_Q(C); 
%     codebook = get_codebook(8); 
%     [v_idx, ~] = get_feedback(Q0, C, codebook); 
%     v = codebook(:, v_idx);
%     v = v/norm(v); 
%     Q0 = Q0 / norm(Q0); 
%     baseline_init = real(v'*Q0'*C*Q0*v/(d1));
%     type_i = type_i + baseline_init; 
%     
%     %% Type II ratio
%     WW = TypeII_rank1(H, H);
%     WW = WW / norm(WW); 
%     type_ii_baseline = real(WW'*C*WW/(d1)); 
%     type_ii = type_ii + type_ii_baseline; 
% end 
% 
% type_i = type_i / n;
% type_ii = type_ii / n;

end