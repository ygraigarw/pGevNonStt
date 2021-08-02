function Nll=GevNll(Prm, X, T);
%function Nll=GevNll(Prm, X, T);
%
%Negative log likelihood for GEV, parameters are linear in time T

Nll=NaN;

%% Linear forms for parameters
Xi=Prm(1)+T*Prm(2);
Sgm=Prm(3)+T*Prm(4);
Mu=Prm(5)+T*Prm(6);

%% Reject non-positive sigma
if min(Sgm)<0;
    return;
end;

%% Reject xi out of sensible range
if min(Xi)<=-0.5 || max(Xi)>0.5; %PhJ20201215
    return;
end;

%% Reject exceedances of upper end point
t0=(1 + Xi.*(X-Mu)./Sgm);
if min(t0)<0;
    return;
end;

%% Calculate NLL
t=t0.^(-1./Xi);
Nll = -sum(-log(Sgm)+(Xi+1).*log(t) - t);

return;