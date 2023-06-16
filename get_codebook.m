function codebook = get_codebook(NP)

% This function is used to get the type_i codebook. 

if NP == 8

    O1 = 4;
    O2 = 4;
    N1 = 2;   % N1 天线的列数
    N2 = 2;   % N2 天线的行数 
    CS_num = N1*N2*2; 
    
    basis = get_basis_(O1, O2, N1, N2, CS_num);
    [~, n] = size(basis);
    
    codebook = [];
    
    for i = 1:n
        v_lm = basis(:, i);
        for j = 1:4
            codeword = [v_lm; exp((2*pi*(j-1)*1i)/4)*v_lm];
            codebook = [codebook, codeword];
        end 
    end 
end 

if NP == 2
    O1 = 4;
    O2 = 1;
    N1 = 2;
    N2 = 1;
    CS_num = N1*N2*2; 
    
    basis = get_basis_(O1, O2, N1, N2, CS_num); 
    codebook = basis; 
end

end

function basis = get_basis_(O1, O2, N1, N2, CS_num)
    for l = 0:O1*N1-1
        for m = 0:O2*N2-1
            u_m = [];
            for i = 0:N2-1
                u_m = [u_m exp(1i*2*pi*m*i/(O2*N2))];
            end
            v_lm = [];
            for j = 0:N1-1
                v_lm = [v_lm exp(1i*2*pi*l*j/(O1*N1))*u_m];
            end
            v_lm = v_lm.';
            basis(:,:,l+1,m+1) = v_lm;
        end  
    end
    basis = reshape(basis, [CS_num/2, O1*O2*N1*N2]);
end