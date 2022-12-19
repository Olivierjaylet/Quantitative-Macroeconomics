% =======================================================================
% LagOrderSelectionVARp.m
% =======================================================================

function results = LagOrderSelectionVARpwilli(y,const,pmax)
% =======================================================================
% Perform and display lag order selection tests for VAR(p) model, i.e.
% Akaike, Schwartz and Hannan-Quinn information criteria
% =======================================================================
% nlag = LagOrderSelectionVARp(y,const,pmax,crit)
% -----------------------------------------------------------------------
% INPUTS
%   - y     : data vector [nobs x nvar]
%   - const : 1 constant, 2 constant+linear trend. [scalar]
%   - pmax  : number of maximum lags to consider. [scalar]
%   - crit  : criteria to compute lag order selection;
%             possible values: 'AIC', 'SIC', 'HQC'
% -----------------------------------------------------------------------
% OUTPUTS
%   - nlag  : number of lags recommended by the selected information crit
% =======================================================================
% Willi Mutschler, December 3, 2021
% willi@mutschler.eu
% =======================================================================

[nobs,K] = size(y); % Sample size and number of variables
T = nobs-pmax;    % Effective sample size used for all estimations, i.e.
                   % number of presample values set aside for estimation
                   % is determined by the maximum order pmax
% Initialize
INFO_CRIT = nan(pmax,3);

% number of presample values set aside for estimation is determined by pmax
% Y = [y_{nlag+1},..., y_{nobs}] is [nvarx(nobs-nlag)] matrix of lagged endogenous variables; note that we need to start in nlag+1 not in 1
Y=transpose(y((pmax+1):nobs,:));
% Z = [Z_{nlag} Z_{nlag+1} ... Z_{nobs-1}] is [(opt.const+nvar*nlag)x(nobs-nlag)] matrix of regressors
Zmax = transpose(lagmatrix(y,[1:pmax]));
Zmax = Zmax(:,pmax+1:nobs); %remove initial observations
% add deterministic terms if any
if const == 1
    Zmax = [ones(1,T); Zmax];
elseif const == 2
    Zmax=[ones(1,T); (pmax+1):nobs; Zmax];
end

for m=1:pmax
    n = const+K*m;                % Number of freely estimated parameters
    Z = Zmax(1:n,:);                 % Data used in estimation
	Ahat = (Y*Z')/(Z*Z');            % OLS and ML estimator
    uhat = Y-Ahat*Z;                 % Residuals
    log_det_SIGMLm = log(det((uhat*uhat')./T)); % ML estimate of variance of errors
    phi_m = (m*K^2+K);
    INFO_CRIT(m,1) = log_det_SIGMLm + 2/T * phi_m;             % Akaike
    INFO_CRIT(m,2) = log_det_SIGMLm + 2*log(log(T))/T * phi_m; % Hannan-Quinn
    INFO_CRIT(m,3) = log_det_SIGMLm + log(T)/T * phi_m;        % Schwartz
end

% Store results and find minimal value of criteria
results = [transpose(1:pmax) INFO_CRIT];