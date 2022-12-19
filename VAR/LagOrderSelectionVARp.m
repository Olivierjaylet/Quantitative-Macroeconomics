% =======================================================================
% LagOrderSelectionARp.m
% =======================================================================

function [nlag, Z_mat, Y] = LagOrderSelectionVARp(ENDO,const,pmin, pmax,crit)
% =======================================================================
% Perform and display lag order selection tests for AR(p) model, i.e.
% Akaike, Schwartz and Hannan-Quinn information criteria
% =======================================================================
% nlag = LagOrderSelectionARp(y,const,pmax,crit)
% -----------------------------------------------------------------------
% INPUTS
%   - ENDO  : data matrix [nobs x nvar]
%   - const : 1 constant, 2 constant+linear trend. [scalar]
%   - pmax  : number of maximum lags to consider. [scalar]
%   - crit  : criteria to compute lag order selection;
%             possible values: 'AIC', 'SIC', 'HQC'
% -----------------------------------------------------------------------
% OUTPUTS
%   - nlag  : number of lags recommended by the selected information crit
% =======================================================================
% Willi Mutschler, November 2, 2021
% willi@mutschler.eu
% =======================================================================

%% construct regressor matrix and dependent variable
% number of presample values set aside for estimation is determined by pmax
%ENDO = ENDO((pmax+1:end),:);
% [T, K]  = size(ENDO);     % sample size

[T, K]  = size(ENDO);     % sample size


nobs_eff = T - pmax;

INFO_CRIT = nan(pmax,1); % initialize
for p=pmin : pmax
    Y = transpose(ENDO((pmax+1):T,:));
    Z = transpose(lagmatrix(ENDO,[1:p]));
    Z = Z(:,pmax+1:T);

    if const == 1
        Z= [ones(1,nobs_eff); Z];
    elseif const == 2
        Z=[ones(1,nobs_eff); (p):nobs_eff+(p-1); Z];
    end
    Z_mat(:,:,p) = Z;
    A = (Y*Z')/(Z*Z'); % OLS and Gaussian ML estimate
    U = Y-A*Z; % OLS and ML residuals
    UUt = U*U'; % sum of squared residuals
    SIGMLu = (1/nobs_eff)*UUt; % ML: not adjusted for # of estimated coefficients

    phi = p * K^2 + K;
    if strcmp(crit,'AIC') % Akaike
        INFO_CRIT(p,:) = log(det(SIGMLu)) + (2/T)*phi;
    elseif strcmp(crit,'SIC') % Schwartz
        INFO_CRIT(p,:) = log(det(SIGMLu)) + (log(T)/T)*phi;
    elseif strcmp(crit,'HQC') % Hannan-Quinn
        INFO_CRIT(p,:) = log(det(SIGMLu)) + (2*log(log(T))/T)*phi;
    end
end

% Store results and find minimal value of criteria
%results = [transpose(1:pmax) INFO_CRIT];
nlag = find(INFO_CRIT == min(INFO_CRIT));
%[~,idx] = sort(results(:,2));
%results = results(idx,:);

% % Display summary of results
% fprintf('*************************************************************\n');
% fprintf('*** OPTIMAL ENDOGENOUS LAGS FROM %s INFORMATION CRITERIA ***\n',crit);
% fprintf('*************************************************************\n');
% disp(array2table(results,'VariableNames',{'Lag',crit}))
% fprintf('  Optimal number of lags (searched up to %d lags):\n',pmax);
% fprintf('  %s Info Criterion:    %d\n',crit,nlag);
% fprintf('\n');