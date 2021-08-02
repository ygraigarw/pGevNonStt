function C=GevMCMC(Y,C);
% function C=GevMCMC(Y,C);
%
% MCMC for GEV fit, all parameters take form p_0 + time * p_1
% Objective is to estimate xi0, xi1, sigma0, sigma1, mu0 and mu1

%% Parameter names
PrmNms={'$\xi_0$';'$\xi_1$';'$\sigma_0$';'$\sigma_1$';'$\mu_0$';'$\mu_1$';};

%% Loop over MCMC iterations
C.AccRat=nan(C.nI,6);
for iI=1:C.nI;
    
    %% Find valid starting solution if iI=1
    if iI==1;
        
        NllStr=GevNll(Y.XSMStart,Y.Dat,Y.Tim);
        
        if isnan(NllStr)==1; %no valid starting solution found. Terminate.
            fprintf(1,'Warning: invalid starting solution.\n');
            fprintf(1,'Trying constant solution.\n');
            tStr=gevfit(Y.Dat);
            Y.XSMStart=[tStr(1);0;tStr(2);0;tStr(3);0];
            NllStr=GevNll(Y.XSMStart,Y.Dat,Y.Tim);
            if isnan(NllStr)==1; %no valid starting solution found. Terminate.
                fprintf(1,'Warning: invalid constant starting solution. Terminating.\n');
                return;
            else;
                Prm=Y.XSMStart;
                Nll=NllStr;
                fprintf(1,'Starting solution found. Starting MCMC\n');
            end;
        else %make the current state the starting state for MCMC
            Prm=Y.XSMStart;
            Nll=NllStr;
            fprintf(1,'Starting solution found. Starting MCMC\n');
        end
        
    end;
    
    %% Iteration counter on screen
    if rem(iI,100)==0;
        fprintf(1,'+');
    elseif rem(iI,10)==0;
        fprintf(1,'.');
    end;
    if rem(iI,1000)==0;
        fprintf(1,'\n');
    end;
    
    %% Loop over parametric forms
    for iP=1:6; %Metropolis Hastings in Gibbs, one parameter at a time
        
        %% Define candidate in terms of current state
        PrmC=Prm;
        
        if iI<=C.AdpItr; %fixed nugget
            PrmC(iP)=PrmC(iP)+randn*C.NgtStr;
        elseif iP==1; %adaptive Metropolis for A B M S
            jP=[1:6]'; nJ=size(jP,1);
            t1=real((1-C.AdpBet)*2.38*sqrtm(cov(C.Prm(max(1,iI-999):iI-1,jP))/nJ)*randn(nJ,1));
            t2=C.AdpBet*0.1*(randn(nJ,1)/nJ);
            PrmC=PrmC+t1+t2;
            C.AccRat(iI-1,jP(2:end))=NaN;
        end;
        
        %% Evaluate likelihood at current state and candidate
        NllC=GevNll(PrmC,Y.Dat,Y.Tim);
        
        if isreal(NllC)==0;
            fprintf(1,'Imaginary NllC\n');
        end
        
        %% MH acceptance step
        if (exp(-NllC+Nll) > rand) && isinf(NllC)==0 && isnan(NllC)==0;
            Prm=PrmC;
            Nll=NllC;
            if iI==1;
                C.AccRat(iI,iP)=1;
            else;
                if iI>100; %Only use last 100 iterations to adjust acceptance rate
                    jI=100;
                else;
                    jI=iI;
                end;
                if iI<=C.AdpItr;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                elseif iP==1;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                else;
                    C.AccRat(iI,iP)=NaN;
                end;
            end;
        else;
            if iI==1; %Only use last 100 iterations to adjust acceptance rate
                C.AccRat(iI,iP)=0;
            else;
                if iI>100; %Only use last 100 iterations to adjust acceptance rate
                    jI=100;
                else;
                    jI=iI;
                end;
                if iI<=C.AdpItr;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                elseif iP==1;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                else;
                    C.AccRat(iI,iP)=NaN;
                end;
            end;
        end;
        
    end;
    
    % Update after complete iteration over variables
    C.Prm(iI,:)=Prm';
    C.Nll(iI)=Nll;
    
    
    if rem(iI,1000)==0;
        
        %% Trace plots
        for j=1:6;
            subplot(2,4,j);
            plot(C.Prm(:,j),'k-');
            title(PrmNms{j},'interpreter','latex');
        end;
        subplot(2,4,7); plot(C.Nll,'k-'); title 'NLL';
        subplot(2,4,8); plot(C.AccRat); title 'Acceptance rates';
        drawnow;
                
    end;
    
    if rem(iI,5000)==0 || iI==C.nI;
        
        %% Save whole chain
        tFil=sprintf('MCMC');
        save(tFil,'C');
        
    end;
    
end;

return;