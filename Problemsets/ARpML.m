function ML = ARpML(y,p,const,alph)

T=size(y,1);

f=@(x) -1*func_LogLikeARpNorm(x,y,p,const);

x0=randn(p+const+1,1);

[x,fval,exitflag,output,grad,hess] = fminunc(f,x0);

thetatilde = x(1:p+const);
sigutilde = x(end);
V = inv(hess);
sd = sqrt(diag(V));
sd_thetatilde = sd(1:p+const);
sd_sigutilde = sd(end);

T_eff = T-p;
logl = -fval;
tstat = thetatilde./sd_thetatilde;
tcrit = -tinv(alph/2,T_eff-p);
pvalues = tpdf(tstat,T_eff-p);

theta_ci=[thetatilde-tcrit.*sd_thetatilde, thetatilde+tcrit.*sd_thetatilde];

ML.T_eff=T_eff;
ML.logl = logl;
ML.thetatilde = thetatilde;
ML.sigutilde = sigutilde;
ML.sd_thetatilde = sd_thetatilde;
ML.sd_sigutilde = sd_sigutilde;
ML.tstat = tstat;
ML.pvalues = pvalues;
ML.theta_ci = theta_ci;

end