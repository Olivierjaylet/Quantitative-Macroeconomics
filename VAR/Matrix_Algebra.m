%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% MATRIX ALGEBRA %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A =[0.5 0 0
    0.1 0.1 0.3
    0 0.2 0.3];
sig_u = [2.25 0 0
    0 1 0.5
    0 0.5 0.74];

%% Eigenvalues of A
EI_1 = eig(A);
% for VAR model, we need to check wether eigen values are lower than one in
% absolute value
% If that is the case for all eigen values, the process is stable and 
% covariance-stationnary
disp(abs(EI_1) <1);
