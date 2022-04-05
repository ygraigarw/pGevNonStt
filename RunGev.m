%% Simple code to fit GEV model, all the parameters of which are assumed to vary linearly in time
% GEV(xi, sgm, mu) with xi=xi_0 + Tim * xi_1, sigma=sigma_0 + Tim * sigma_1, mu=mu_0 + Tim * mu_1, 
% Will run as is for toy data to check 
% Need to input data as structure (see occurrences of USER INPUT) below

%% Set up
clc; clear; clf; pLtx;
VrbNms={'$\xi$';'$\sigma$';'$\mu$'};

%% Simulate a sample of data
if 1; %for testing

    % Time variable (assumed to be defined on 0,1)
    X.nT=10000;
    X.Tim=linspace(0,1,X.nT)';
    
    % Block length in years (e.g. X.BlcLngYrs=5 means that we are looking at 5-year maxima data)
    % If BlcLngYrs is not defined, assume it is =1
    X.BlcLngYrs=5;
    
    % True parameters P0=[xi0;xi1;sgm0;sgm1;mu0;mu1] of linear regression 
    if 1; %Obvious difference
    X.Prm0=[...
        -0.3;0.1;
        1;1;...
        5;1;...
        ];
    else; %no difference
    X.Prm0=[...
        -0.3;0;
        1;0;...
        5;0;...
        ];
    end;
    % e.g. mu0 is the location parameter at T=0 (start), and mu0+mu1 is the location parameter at T=1 (end)
    
    % True parameter estimates in time
    X.XSM0=[ones(X.nT,1)*X.Prm0(1)+X.Tim*X.Prm0(2) ones(X.nT,1)*X.Prm0(3)+X.Tim*X.Prm0(4) ones(X.nT,1)*X.Prm0(5)+X.Tim*X.Prm0(6)];
    
    % Generate data from GEV
    X.Dat=gevrnd(X.XSM0(:,1),X.XSM0(:,2),X.XSM0(:,3));
    
    X, % See the structure
    
end;

%% ***USER INPUT*** Read in your data here
if 0; %set to zero if you want to use simulated data from above
     %X.nT ;  %  1 x 1 number of time points
     %X.Tim ; % nT x 1 times on [0,1]
     %X.Dat ; % nT x 1 data 
     load userInput.mat;
end;

%% Find starting solution by regression
if 1; 
    
    % Split the data into blocks
    % Fit a stationary GEV model per block
    % Fit straight lines through the parameter estimates per block using linear regression
    
    Y.nT=X.nT;
    Y.Tim=X.Tim;
    Y.Dat=X.Dat;
    Y.nT=size(X.Dat,1);
    Y.nB=10;
    Y.Blc=pMakCV(Y.nT,Y.nB,Y.Tim);
    
    % Check for existence of X.BlcLngYrs; if it exists, use it; otherwise set Y.BlcLngYrs=1
    if isfield(X,'BlcLngYrs')==1;
        Y.BlcLngYrs=X.BlcLngYrs;
    else;
        Y.BlcLngYrs=1;
    end;
    
    clf;
    subplot(2,1,1); plot(Y.Tim, Y.Dat,'k.');
    xlabel 'Time';
    ylabel 'Value';
    pAxsLmt;pDflBig;
    
    tRgr=nan(Y.nB,4);
    for iB=1:Y.nB;
        tPrm=gevfit(Y.Dat(Y.Blc==iB));
        tTim=mean(Y.Tim(Y.Blc==iB));
        tRgr(iB,:)=[tTim tPrm];
    end;
    [jnk,tOrd]=sort(tRgr(:,1));
    Y.Rgr=tRgr(tOrd,:);
    
    Y.XSMStart=nan(6,1);
    for j=1:3;
        subplot(2,3,3+j); plot(Y.Rgr(:,1),Y.Rgr(:,j+1),'b.-');
        tPrm=[ones(Y.nB,1) Y.Rgr(:,1)]\Y.Rgr(:,j+1);
        Y.XSMStart(2*j-1:2*j)=tPrm';
        subplot(2,3,3+j); hold on; plot(Y.Rgr(:,1),ones(Y.nB)*Y.XSMStart(2*j-1)+Y.Rgr(:,1)*Y.XSMStart(2*j),'k.-');
        ylabel(VrbNms{j},'interpreter','latex'); 
        if j==1; xlabel 'Time'; end;
        pAxsLmt;pDflBig;
    end;
    
end;

%% ***USER INPUT*** Find better solution using MCMC
if 1;
    
    C.nI=10000;       % Number of MCMC iterations - 1e4 minimum when used in anger
    C.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
    
    C.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
    C.AdpItr=1000;    % Number of warm up iterations - don't change
    C.AdpBet=0.05;    % Adaptive MC - don't change
    
    C=GevMCMC(Y,C);   % Run MCMC algorithm
    
end;

%% Plot MCMC results
if 1;
    
    clf;
    PrmNms={'$\xi_0$';'$\xi_1$';'$\sigma_0$';'$\sigma_1$';'$\mu_0$';'$\mu_1$';};
    load MCMC;
    for j=1:6;
        subplot(2,4,j);
        pHst(C.Prm(C.nI-C.n2Plt+1:end,j));
        title(PrmNms{j},'interpreter','latex');
        pAxsLmt; pDflBig;
        fprintf(1,'Median %s = %g\n',PrmNms{j},quantile(C.Prm(C.nI-C.n2Plt+1:end,j),0.5));
    end;
    subplot(2,4,7); plot(C.Nll,'k-'); title 'NLL'; pAxsLmt; pDflBig;
    subplot(2,4,8); plot(C.AccRat); title 'Acceptance rates'; pAxsLmt; pDflBig;
    pGI('GevNonStt-McmcTracePlot',2);
    
end;


%% Compare 100-year values at start and end of study as a function of threshold
% Added 20201201 for consistency with GP code
if 1;

    %% Calculate return values
    C.RV.RtrPrd=100; % Return period of interest
    C.RV.nRls=1000;  % Number of realisations to use (1000 is good)
    
    % Parameter estimates at start
    t=randi(C.nI-C.n2Plt,C.RV.nRls,1)+C.n2Plt;
    tXi=C.Prm(t,1);
    tSgm=C.Prm(t,3);
    tMu=C.Prm(t,5);
    C.RV.Est(:,1)=(tSgm./tXi).*( (-log(1-Y.BlcLngYrs/C.RV.RtrPrd)).^(-tXi) - 1 ) + tMu; % Return value for year zero
        
    % Parameter estimates at end
    tXi=C.Prm(t,1)+C.Prm(t,2);
    tSgm=C.Prm(t,3)+C.Prm(t,4);
    tMu=C.Prm(t,5)+C.Prm(t,6);
    C.RV.Est(:,2)=(tSgm./tXi).*( (-log(1-Y.BlcLngYrs/C.RV.RtrPrd)).^(-tXi) - 1 ) + tMu; % Return value for year zero
    
    % Summary statistics of differences
    C.RV.Prb=mean(C.RV.Est(:,2)>C.RV.Est(:,1)); % Prob. that RVEnd > RVStart
    C.RV.Cdf=[sort(C.RV.Est(:,1)) sort(C.RV.Est(:,2))];
    C.RV.Qnt=[quantile(C.RV.Est(:,1),[0.025 0.5 0.975]) quantile(C.RV.Est(:,2),[0.025 0.5 0.975])]; % Quantiles
    
    
    %% Figure
    clf; 
    C.RV.PrbVls=((1:C.RV.nRls)'-0.5)/C.RV.nRls;
    subplot(1,2,1); hold on;
    plot(C.RV.Cdf(:,1),C.RV.PrbVls,'k');
    plot(C.RV.Cdf(:,2),C.RV.PrbVls,'r');
    if isfield(X,'Prm0')==1; % True values are known
        tXi=X.Prm0(1);
        tSgm=X.Prm0(3);
        tMu=X.Prm0(5);
        tRVTru=(tSgm./tXi).*( (-log(1-Y.BlcLngYrs/C.RV.RtrPrd)).^(-tXi) - 1 ) + tMu;
        plot(tRVTru*ones(2,1),[0 1]','k--','linewidth',2);
        tXi=X.Prm0(1)+X.Prm0(2);
        tSgm=X.Prm0(3)+X.Prm0(4);
        tMu=X.Prm0(5)+X.Prm0(6);
        tRVTru=(tSgm./tXi).*( (-log(1-Y.BlcLngYrs/C.RV.RtrPrd)).^(-tXi) - 1 ) + tMu;
        plot(tRVTru*ones(2,1),[0 1]','r--','linewidth',2);
    end;
    pAxsLmt; pDflHug;
    title 'Distribution of RVStart (k), RVEnd (r) [True - - -]'; 
    xlabel(sprintf('%g-year maximum value',Y.BlcLngYrs));
    ylabel(sprintf('$F_{{%gYearMaximum}}$',Y.BlcLngYrs),'interpreter','latex');
    pAxsLmt; pDflHug;
    %
    subplot(1,2,2); hold on;
    plot(C.RV.Cdf(:,1),log10(1-C.RV.PrbVls),'k','linewidth',2);
    plot(C.RV.Cdf(:,2),log10(1-C.RV.PrbVls),'r','linewidth',2);
    pAxsLmt; pDflHug;
    title 'Distribution of RVStart (k), RVEnd (r) [True - - -]'; 
    xlabel(sprintf('%g-year maximum value',Y.BlcLngYrs));
    ylabel(sprintf('$\\log_{10}(1-F_{{%gYearMaximum}})$',Y.BlcLngYrs),'interpreter','latex');
    pAxsLmt; pDflHug;
    %
    pDatStm(sprintf('Prb(RVEnd$>$RVStart)=%4.2f\n',C.RV.Prb)); pGI('RV',2);
    pGI('GevNonStt-Compare100YearReturnValues',2);
    
    %% Summary statistics to screen
    clc;
    fprintf(1,'SUMMARY\n');
    fprintf(1,'Quantiles 2.5pc 50pc 97.5pc START: %g %g %g\n',C.RV.Qnt(1:3));
    fprintf(1,'Quantiles 2.5pc 50pc 97.5pc END: %g %g %g\n',C.RV.Qnt(4:6));
    if C.RV.Prb>0.975 || C.RV.Prb<=0.025;
        fprintf(1,'Prb(RVEnd>RVStart)=%4.2f SIGNIFICANT\n',C.RV.Prb);
    else;
        fprintf(1,'Prb(RVEnd>RVStart)=%4.2f (not significant)\n',C.RV.Prb);
    end;

end;