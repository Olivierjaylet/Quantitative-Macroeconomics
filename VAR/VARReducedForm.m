%%%% OLS Function for VAR model %%%%

%%% Input %%%
% ENDO : [KxT] data Matrix [nobs x nvar]
% nlag : [integer] number of lags
% const : [flag] : 0 no constant; 1 constant; 2 constant + linear trend


%%% Output %%%
% VAR : Structure :
%   - T_eff :       [scalar]        effective sample size used in estimation
%   - thetathat :   [(const+p)x1]   estimation of coefficients
%   - siguhat :     [scalar]        estimate of standard deviation of error term
%   - sd_thetahat : [(const+p)x1]   estimate of standard error of coefficients
%   - tstats :      [(const+p)x1]   t statistics
%   - pvalues :     [(const+p)x1]   p values of H_0 : thetahat=0
%   - theta_ci :    [(const+p)x2]   (1-alph)% confidence intervall for theta given
%                                   significance level alph
%   - resid :       [T_effx1]       residuals



function VAR = VARReducedForm(ENDO,nlag,const)
[nobs, nvar]=size(ENDO);                % sample size
nobs_eff = nobs-nlag;                   % effective sample size used in estimation

% Create Y matrix (need transpose of ENDO and start at nlag + 1)
Y = transpose(ENDO((nlag+1):nobs,:));

% Create Z matrix
Z = transpose(lagmatrix(ENDO, [1 : nlag]));
Z = Z(:,nlag+1:nobs); % Remove initial observations

if const ==1
    Z = [ones(1, nobs_eff); Z];             % constant term
elseif const == 2
    Z = [ones(1, nobs_eff); (nlag+1):obs; Z];    % constant term + time trend
end



A = (Y*Z')/(Z*Z');   % OLS estimate;
U =  Y-Y*A;          % residuals
UUt = U*U';          % SSR 

SIGOLSu = (1/(nobs_eff-nvar*nlag-opt.const))*UUt; % OLS: adjusted for # of estimated coefficients
SIGMLu = (1/nobs_eff)*UUt; % ML: not adjusted for # of estimated coefficients
% Compute maximum absolute Eigenvalue of companion VAR(1) A matrix to check for stability

Acomp = [A(:,1+opt.const:nvar*nlag+opt.const);
eye(nvar*(nlag-1)) zeros(nvar*(nlag-1),nvar)];
maxEig = max(abs(eig(Acomp)));

%% Save into structure
VAR.ENDO = ENDO;
VAR.nlag = nlag;
VAR.opt = opt;
VAR.Z = Z;
VAR.Y = Y;
VAR.A = A;
VAR.residuals = U;
VAR.SigmaOLS = SIGOLSu;
VAR.SigmaML = SIGMLu; % Maximum Likelihood COV Matrix is not adjusted for # of estimated coefficients
VAR.Acomp = Acomp;
VAR.maxEig = maxEig;

end