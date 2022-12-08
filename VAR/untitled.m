y=importdata('C:\Users\User\OneDrive\Bureau\Quantitative-Macroeconomics\VAR\USOil.csv');
y=y.data;
nlag = 4;
opt.const = 1;
VAR = VARReducedForm(y,nlag,opt);
SIGu = VAR.SigmaOLS;
B0inv = chol(SIGu,'lower'); % Using choleski
disp(inv(B0inv));
%%
opt.nsteps = 20;
IRFs(VAR.Acomp, BOinv, opt);