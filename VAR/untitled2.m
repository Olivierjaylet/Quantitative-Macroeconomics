maxLag = 12;
ENDO=ENDO(1+maxLag:end,:);
[T,K] = size(ENDO);

for nlag=1:maxLag
    nobs_eff = nobs - nlag; % Effective sample size used in estimation
    
    %% Create independent vector and lagged dependent matrix
    % Y = [y_{nlag+1},..., y_{nobs}] is [nvarx(nobs-nlag)] matrix of lagged endogenous variables; note that we need to start in nlag+1 not in 1
    Y = transpose(ENDO((nlag+1):nobs,:));
    
    % Z = [Z_{nlag} Z_{nlag+1} ... Z_{nobs-1}] is [(opt.const+nvar*nlag)x(nobs-nlag)] matrix of regressors
    Z = transpose(lagmatrix(ENDO,[1:nlag]));
    Z = Z(:,lagmax+1:nobs); %remove initial observations
    % add deterministic terms if any
    if opt.const == 1
        Z = [ones(1,nobs_eff); Z];
    elseif opt.const == 2
        Z=[ones(1,nobs_eff); (nlag+1):nobs; Z];
    end
    
    %% Compute the matrix of coefficients & Covariance Matrix
    A = (Y*Z')/(Z*Z'); % OLS and Gaussian ML estimate
    U = Y-A*Z; % OLS and ML residuals
    UUt = U*U'; % sum of squared residuals
    SIGMLu = (1/nobs_eff)*UUt; % ML: not adjusted for # of estimated coefficients
    phi = p * K^2 + K;
    if strcmp(crit,'AIC') % Akaike
        INFO_CRIT(p,:) = log(det(SIGMLubar)) + (2/T)*phi;
    elseif strcmp(crit,'SIC') % Schwartz
        INFO_CRIT(p,:) = log(det(SIGMLubar)) + (log(T)/T)*phi;
    elseif strcmp(crit,'HQC') % Hannan-Quinn
        INFO_CRIT(p,:) = log(det(SIGMLubar)) + (2*log(log(T))/T)*phi;
    end
end