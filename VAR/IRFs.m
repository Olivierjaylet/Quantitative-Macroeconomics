function IRFpoint  = IRFs(Acomp,impact,nsteps,IRFcumsum,varnames,epsnames,noplot)
% =======================================================================
% Compute IRFs for a SVAR model
% B0*y_t = B_1*y(-1) + B_2*y(-2) +...+ B_p*y(-p) + e_t with e_t ~ iid(0,I)
% given companion matrix A of the reduced-form VAR and impact matrix inv(B0)
% =======================================================================
% IRFpoint  = IRFs(Acomp,impact,nsteps,IRFcumsum,vnames,epsnames,doplot,dosave,filename)
% -----------------------------------------------------------------------
% INPUT
%   - Acomp:       [(nvars*(nlags-1)) x (nvars*nlags)]  companion matrix of reduced-form VAR
%   - impact:      [nvars x nvars]                      impact matrix of SVAR model, either inv(B0) or inv(B0)*sqrt(E[epsilon_t*epsilon_t'])
%   - IRFcumsum:   [nvars x 1]                          boolean vector, indicating for which variable to compute cumulative sum in IRFs (1)
%   - varnames:    [nvars x 1]                          vector of strings with variable names
%   - epsnames:    [nvars x 1]                          vector of strings with structural shock names
%   - noplot       boolean                              1: turn off displaying of plots (useful for bootstrapping)
% ----------------------------------------------------------------------- 
% OUTPUT
%   - IRFpoint(j,k,h)  [nvars,nvars,nsteps+1]         array with h=0,1,...,nsteps+1 steps, containing the IRF of the 'j' variable to the 'k' shock
% =======================================================================
% Willi Mutschler, December 8, 2022
% willi@mutschler.eu

nvars = size(impact,1);
nlag  = size(Acomp,2)/nvars;

% set default options if not specified
if nargin < 4 || isempty(IRFcumsum)
    IRFcumsum = false(1,nvars);
end
if nargin < 5 || isempty(varnames)
    varnames = "y" + string(1:nvars);
end
if nargin < 6 || isempty(epsnames)
    epsnames = "e" + string(1:nvars);
end
if nargin < 7 || isempty(noplot)
    noplot = false;
end

% initialize variables
IRFpoint = nan(nvars,nvars,nsteps+1);
J = [eye(nvars) zeros(nvars,nvars*(nlag-1))];

% compute the impulse response function
Ah = eye(size(Acomp)); % A^0
JtB0inv = J'*impact;
for h=1:(nsteps+1)
    IRFpoint(:,:,h) = J*Ah*JtB0inv;
    Ah = Ah*Acomp; %A^h = A^(h-1)*A
end

% use cumsum to get response of level variables from original variables in differences
for ivar = 1:nvars
    if IRFcumsum(ivar) == true
        IRFpoint(ivar,:,:) = cumsum(IRFpoint(ivar,:,:),3);
    end
end

% plot
if ~noplot
    % Define a timeline
    steps = 0:1:nsteps;
    x_axis = zeros(1,nsteps+1);
    figure('units','normalized','outerposition',[0 0.1 1 0.9]);
    count = 1;
    for ivars=1:nvars % index for variables
        for ishocks=1:nvars % index for shocks    
            irfpoint = squeeze(IRFpoint(ivars,ishocks,:));
            subplot(nvars,nvars,count);
            plot(steps,irfpoint,steps,x_axis,'k','LineWidth',2);
            xlim([0 nsteps]);
            title(epsnames(ishocks), 'FontWeight','bold','FontSize',10);
            ylabel(varnames(ivars), 'FontWeight','bold','FontSize',10);
            count = count+1;
        end
    end
end