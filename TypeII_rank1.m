
function [WW] = TypeII_rank1(WB_H,SB_H)
Nr = 2;
Nt = 32;

N_PSK = 8;
rank = 1;
O1 = 4;
O2 = 4;
N1 = 8;
N2 = 2;
P_CSI_RS = 32;
L = 4;

% -----------------------------------------
% Quantization Level
subband_phase_list = (0:8)*(2*pi)/N_PSK;

wideband_amplitude_list = [0];
for i = 6:-1:0
    wideband_amplitude_list = [wideband_amplitude_list sqrt(1/2^i)];
end

subband_amplitude_list = [sqrt(1/2) 1];

% -----------------------------------------
% Get the O1*O2(16) orthogonal_basis with N1*N2 (N1*N2,1) vector
num = 1;
for q1 = 0:O1-1
    for q2 = 0:O2-1      
        orthogonal_basis = [];
        for n1 = 0:N1-1
            l = O1*n1+q1;
            for n2 = 0:N2-1
                m = O2*n2+q2;
                u_m = [];
                for i = 0:N2-1
                    u_m = [u_m exp(1i*2*pi*m*i/(O2*N2))];
                end
                v_lm = [];
                for j = 0:N1-1
                    v_lm = [v_lm exp(1i*2*pi*l*j/(O1*N1))*u_m];
                end
                v_lm = v_lm.';
                orthogonal_basis = [orthogonal_basis v_lm];
            end
        end
        orthogonal_basis_set(:,:,num) = orthogonal_basis;
        num = num+1;
    end
end


% Add up each correlation matrix to get WB (600 subcarrier(50 RB))

[V_eig,D_eig] = eig(WB_H'* WB_H/100);
[U,S,V] = svd(WB_H);
D_eig = diag(D_eig);
[~,I_eig]=sort(abs(D_eig),'descend');

% For rank=2, Primary and secondary eigenvectors
V_eig = V_eig(:,I_eig(1:rank));


% For WB, determine the orthogonal basis
projection_value_list = [];

for index = 1:O1*O2
    orthogonal_basis = orthogonal_basis_set(:,:,index);
    projection_list = [];
    for k = 1:N1*N2
        projection = 0;
        for ll=1:rank
            d1 = V_eig(1:P_CSI_RS/2,ll);
            d2 = V_eig(P_CSI_RS/2+1:end,ll);
            projection1 = d1'*orthogonal_basis(:,k)/norm(orthogonal_basis(:,k));
            projection2 = d2'*orthogonal_basis(:,k)/norm(orthogonal_basis(:,k));
            projection = projection + abs(projection1)^2+abs(projection2)^2;
        end
        projection_list = [projection_list projection];
    end
    [value,projection_index] = sort(projection_list,'descend');
    projection_value = sum(value(1:L));
    projection_value_list = [projection_value_list projection_value];
    candidate_orthogonal_basis_vector(:,:,index) = orthogonal_basis(:,projection_index(1:L));
end

[~, max_index] = max(projection_value_list);
orthogonal_basis_vector = candidate_orthogonal_basis_vector(:,:,max_index);



% ll denotes rank
for ll=1:rank
    d1 = V_eig(1:P_CSI_RS/2,ll);
    d2 = V_eig(P_CSI_RS/2+1:end,ll);
     
    % Get WB amplitude
    for i = 1:L
        projection_i1 = orthogonal_basis_vector(:,i)' * d1/norm(orthogonal_basis_vector(:,i));
        p_l_1(i,ll) = abs(projection_i1);

        projection_i2 = orthogonal_basis_vector(:,i)' * d2/norm(orthogonal_basis_vector(:,i));
        p_l_1(i+L,ll) = abs(projection_i2);
    end
    
    % Normalization WB amplitude
    max_p_l_1(ll) = max(p_l_1(:,ll));
    p_l_1(:,ll) = p_l_1(:,ll)/max_p_l_1(ll);
    
    % Quantization WB amplitude
    for i = 1:L
        [~, ii] = min(abs(p_l_1(i,ll)-wideband_amplitude_list'));
        p_l_1(i,ll) = wideband_amplitude_list(ii);
        
        [~, ii] = min(abs(p_l_1(i+L,ll)-wideband_amplitude_list'));
        p_l_1(i+L,ll) = wideband_amplitude_list(ii);      
    end
end

% Add up 4 correlation matrix to get SB and calculate eigenvectors
[V_eig,D_eig] = eig(SB_H'* SB_H/4);
D_eig = diag(D_eig);
[~,I_eig]=sort(abs(D_eig),'descend');

V_eig = V_eig(:,I_eig(1:2));

% --------------------------------------
% Get SB amplitude and phase
for ll = 1:rank
    d1 = V_eig(1:P_CSI_RS/2,ll);
    d2 = V_eig(P_CSI_RS/2+1:end,ll);
    
    for i = 1:L
        projection_i1 = orthogonal_basis_vector(:,i)' * d1/norm(orthogonal_basis_vector(:,i));
        phase = angle(projection_i1);
        phase = mod(phase,2*pi);
        
        [~, ii] = min(abs(phase-subband_phase_list')); 
        phi_l(i,ll) = subband_phase_list(ii);
            
        amplitude = abs(projection_i1);
        amplitude = amplitude/max_p_l_1(ll);
        amplitude = amplitude/p_l_1(i,ll);
        [~, ii] = min(abs(amplitude-subband_amplitude_list'));
        p_l_2(i,ll) = subband_amplitude_list(ii);
        
%         ------------------
        projection_i2 = orthogonal_basis_vector(:,i)' * d2/norm(orthogonal_basis_vector(:,i));
        phase = angle(projection_i2);
        phase = mod(phase,2*pi);
        
        [~, ii] = min(abs(phase-subband_phase_list')); 
        phi_l(i+L,ll) = subband_phase_list(ii);
            
        amplitude = abs(projection_i2);
        amplitude = amplitude/max_p_l_1(ll);
        amplitude = amplitude/p_l_1(i+L,ll);
        [~, ii] = min(abs(amplitude-subband_amplitude_list'));
        p_l_2(i+L,ll) = subband_amplitude_list(ii);
   
    end
end

% --------------------------------------
% Calculate W
for ll = 1:rank
    parameter = 0;
    for i = 1:2*L
        parameter = parameter + (p_l_1(i,ll) * p_l_2(i,ll))^2;
    end
    parameter = sqrt(N1*N2*parameter);
    
    up_part = zeros(P_CSI_RS/2,1);
    down_part = zeros(P_CSI_RS/2,1);
    for i = 1:L
        up_part = up_part + p_l_1(i,ll) * p_l_2(i,ll) * exp(1j * phi_l(i,ll)) * orthogonal_basis_vector(:,i);
        down_part = down_part + p_l_1(i+L,ll) * p_l_2(i+L,ll) * exp(1j * phi_l(i+L,ll)) * orthogonal_basis_vector(:,i);
    end
    
    vector  = [up_part; down_part];
    
    W(:,ll) = parameter * vector;

end

% WW = 1/sqrt(2)*[W(:,1) W(:,2)];
WW = W;
WW = WW/norm(WW) * norm(SB_H);

end