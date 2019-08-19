function [V0,Fd0,Fdexpected,FdsumP,Optionval,Optdecision,Projectorder,FixStoprule,runtime]=ANDPfinal(t0fix,A,q,Qtp,Yset,Yprime,Linkprime,D,lambda,ro1,ro2,T,P,numBasis)

%t0fix is Nx1 vector of initial travel times on each link
%A is MxN node-link incidence matrix
%q is MxM matrix of OD trip table for the given period
%D is Nx1 matrix for the construction cost per link
%B is a scalar value for budget
%Yset is Nx1 vector of committed link investment design
%Linkprime is an N'xH matrix of H cols of projects and N' link investments
%T is a scalar time horizon for planning
%P is the scalar number of scenarios to simulate
%NOSP originated from RONDP9b - this computes the option value of staging
%links - order fixed - NOSP1 will compute the max of these
%12/12/10 - ro1 is risk free, ro2 is risk adjusted

tic;
numlinks=size(t0fix,1);
numlinkinv=size(Yprime,1);
numprojects=size(Linkprime,2);
savedFdsumP=cell(T+1,numprojects);
indexFdsumP=cell(T+1,numprojects);

Dfix=zeros(numlinks,1);    %creating a D of size(A,2)
for i=1:numlinkinv
    Dfix(Yprime(i,1),1)=D(i,1);
end;

Xout0d=zeros(numlinks,T,P);

for p=1:P
    for t=1:T       %generate the flows
        Xout0d(:,t,p)=AON7(t0fix,A,Qtp(:,:,t,p),[]);    
    end
end
Xout00=AON7(t0fix,A,q,[]);

%initiate with 1st project
Projectorder=1;     %defines current list of projects as the 1st project

for k=2:numprojects
    Bestval=-inf;
    for r=1:k     %this is for each possible ordering of additional project to the current list (Projectorder)
        %%need to update extra row
        Optionval=zeros(k,1);
        Optdecision=int8(zeros(k,1));
        FdsumP=zeros(P,T,k);
        Fd0=zeros(k,1);
        V0=zeros(k,1)+realmin;
        Fdexpected=zeros(k,1);  %for the expected deferral option at time 0
        FixStoprule=int8(zeros(T,P,k));        
        if r==1
            Temporder=[k Projectorder];
            Bestorder=Temporder;
        elseif r==k
            Temporder=[Projectorder k];
        else
            Temporder=[Projectorder(1:r-1) k Projectorder(r:k-1)];
        end
        Yfix0=zeros(numlinks,k);
        Fdsum0=zeros(k,1);
        Linkfix=zeros(numlinks,k);
        for i=1:k
            s=1;
            stop=0;
            while and(s<=size(indexFdsumP{1,i},1),stop==0)
                if sort(indexFdsumP{1,i}(s,:))==sort(Temporder(1,1:i))
                    Fdsum0(i,1)=savedFdsumP{1,i}(s,1);
                    stop=1;
                end
                s=s+1;
            end
            if Fdsum0(i,1)==0
                Ysetprime=zeros(numlinkinv,1);
                Linksetprime=zeros(numlinkinv,1);
                for u=1:i
                    Ysetprime=Ysetprime+Yset.*Linkprime(:,Temporder(1,u));
                    Linksetprime=Linksetprime+Linkprime(:,Temporder(1,u));
                end
                for j=1:numlinkinv
                    Yfix0(Yprime(j,1),i)=Ysetprime(j,1);
                    Linkfix(Yprime(j,1),i)=Linksetprime(j,1);
                end;
                timeY=t0fix.*(1-Yfix0(:,i));
                XoutY=AON7(timeY,A,q,[]); %this calls an All-or-Nothing assignment, it can be replaced by any other type of traffic assignment procedure (UE, SO, etc)
                Vtemp=(Xout00'*t0fix-XoutY'*timeY)*lambda;
                Vtemp=Vtemp/ro2;   %fixed this in GOLIDOSq but not in any of the NIDO/NODP models, 6/13/10
                V0(i,1)=Vtemp;  
                Fdsum0(i,1)=V0(i,1)-Linkfix(:,i)'*Dfix;
                indexFdsumP{1,i}=[indexFdsumP{1,i}; Temporder(1,1:i)];
                savedFdsumP{1,i}=[savedFdsumP{1,i}; Fdsum0(i,1)];
            end
        end

        %%this is where we we would break out into link combinations of
        %%numprojects

        Fdmean=zeros(P,k);
        Vd=zeros(T,P,k);  %changing this to the sum of links instead of link-benefits
        t=T;
        while t>0     %calculate the Ft's
            for i=1:k
                s=1;
                stop=0;
                while and(s<=size(indexFdsumP{1+t,i},1),stop==0)
                    if sort(indexFdsumP{1+t,i}(s,:))==sort(Temporder(1,1:i))
                        FdsumP(:,t,i)=savedFdsumP{1+t,i}(s,:)';
                        stop=1;
                    end
                    s=s+1;
                end
                if FdsumP(:,t,i)==0
                    for p=1:P
                        timeYd=t0fix.*(1-Yfix0(:,i));
                        XoutYd=AON7(timeYd,A,Qtp(:,:,t,p),[]);
                        Vdtemp=(Xout0d(:,t,p)'*t0fix-XoutYd'*timeYd)*lambda;
                        Vdtemp=Vdtemp/ro2;
                        Vd(t,p,i)=Vdtemp;
                        FdsumP(p,t,i)=Vd(t,p,i)-Linkfix(:,i)'*Dfix;
                    end
                    indexFdsumP{t+1,i}=[indexFdsumP{t+1,i}; Temporder(1,1:i)];
                    savedFdsumP{t+1,i}=[savedFdsumP{t+1,i}; FdsumP(:,t,i)'];
                end
            end
            t=t-1;
        end

        %this is for M = 5
        t=T;
        Fd=zeros(P,T,k);  %individual link + option value
        RegScale=mean(V0(2:k,1)-V0(1:k-1,1));
        while t>0
            for h=0:k-1
                dbasis=zeros(P,2); 
                RegYd=zeros(P,1);  
                kd=0; 
                for p=1:P
                    if t==T
                        if h==0
                            Fd(p,t,k-h)=max(0,FdsumP(p,t,k-h)-FdsumP(p,t,k-h-1));
                        elseif k-h==1
                            Fd(p,t,k-h)=max(0,FdsumP(p,t,k-h)+Fd(p,t,k-h+1));
                        else
                            Fd(p,t,k-h)=max(0,FdsumP(p,t,k-h)-FdsumP(p,t,k-h-1)+Fd(p,t,k-h+1));
                        end
                        if Fd(p,t,k-h)>0
                            FixStoprule(T,p,k-h)=1;
                        else
                            FixStoprule(T,p,k-h+1:k)=0;                        
                        end
                    else       
                        if h==0
                            Fd(p,t,k-h)=FdsumP(p,t,k-h)-FdsumP(p,t,k-h-1);
                        elseif k-h==1
                            Fd(p,t,k-h)=FdsumP(p,t,k-h)+Fd(p,t,k-h+1);
                        else
                            Fd(p,t,k-h)=FdsumP(p,t,k-h)-FdsumP(p,t,k-h-1)+Fd(p,t,k-h+1);
                        end
                        if Fd(p,t,k-h)>0
                            kd=kd+1;
                            dbasis(kd,1)=Fd(p,t,k-h);  %first column is the value to regress
                            dbasis(kd,2)=p;       %second column is the path to update later
                            RegYd(kd,1)=Fd(p,t+1,k-h)/(1+ro2)+Linkfix(:,k-h)'*Dfix*(1/(1+ro2)-1/(1+ro1));
                        else
                            Fd(p,t,k-h)=0;   
                            FixStoprule(t,p,k-h+1:k)=0;
                        end
                    end

                end
                RegYd=RegYd(1:kd);
                dbasis=dbasis(1:kd,:);
                RegYd=RegYd/RegScale;
                dbasis(:,1)=dbasis(:,1)/RegScale;

                if t<T
                    if kd>0
                        RegXd=zeros(kd,numBasis+1);
                    end
                    for j=1:kd   %recursive generation of Hermite polynomial basis functions
                        RegXd(j,1)=1;
                        RegXd(j,2)=dbasis(j,1);
                        for m=2:numBasis
                            RegXd(j,m+1)=dbasis(j,1)*RegXd(j,m)-m*RegXd(j,m-1);
                        end
                    end
                    if kd>0
                        dCoeff=((RegXd'*RegXd)\RegXd')*RegYd;   %least squares regression
                    end
                    for j=1:kd
                        if isnan(dCoeff(1,1))
                            if Fd(dbasis(j,2),t,k-h)>mean(RegYd)
                                FixStoprule(t,dbasis(j,2),k-h)=1;
                            else
                                Fd(dbasis(j,2),t,k-h)=mean(RegYd);
                                FixStoprule(t,dbasis(j,2),k-h+1:k)=0;
                            end                            
                        else
                            if Fd(dbasis(j,2),t,k-h)>RegScale*RegXd(j,:)*dCoeff
                                FixStoprule(t,dbasis(j,2),k-h)=1;
                            else
                                Fd(dbasis(j,2),t,k-h)=RegScale*RegXd(j,:)*dCoeff;
                                FixStoprule(t,dbasis(j,2),k-h+1:k)=0;
                            end
                        end
                    end

                end
            end
            t=t-1;
        end

        for h=1:k
            for p=1:P
                t=1;
                while t<=T
                    if FixStoprule(t,p,h)==1
                        Fdmean(p,h)=Fd(p,t,h)/(1+ro2)^t+Linkfix(:,h)'*Dfix*(1/(1+ro2)^t-1/(1+ro1)^t);
                        t=T+1;
                    else
                        t=t+1;
                    end
                end
            end
        end

        for h=0:k-1
            Fdexpected(k-h,1)=mean(Fdmean(:,k-h));
            if k-h==1
                Fd0(k-h,1)=Fdsum0(k-h,1);
            else
                Fd0(k-h,1)=Fdsum0(k-h,1)-Fdsum0(k-h-1,1);
            end
            if h==0
                if Fd0(k-h,1)>Fdexpected(k-h,1)
                    Optdecision(k-h,1)=1;
                    Optionval(k-h,1)=Fd0(k-h,1);
                else
                    Optionval(k-h,1)=Fdexpected(k-h,1);
                end
            else
                if Fd0(k-h,1)+Optionval(k-h+1,1)>Fdexpected(k-h,1)
                    Optdecision(k-h,1)=1;
                    Optionval(k-h,1)=Fd0(k-h,1)+Optionval(k-h+1,1);
                else
                    Optionval(k-h,1)=Fdexpected(k-h,1);
                    Optdecision(k-h+1:k,1)=0;
                end
            end
        end
        %the following chooses the best insertion point for the new project
        Tempval=Optionval(1,1);  %11/16/10 - fixed to new criteria
        if Tempval>Bestval
            Bestval=Tempval;
            Bestorder=Temporder;
            BestOptionval=Optionval;
            Bestdecision=Optdecision;
            BestFdsumP=FdsumP;
            BestFd0=Fd0;
            BestV0=V0;
            BestFdexpected=Fdexpected;
            BestFixStoprule=FixStoprule;
        end
        
    end
    Projectorder=Bestorder;
    
    if k==numprojects
        %at the end we assign the values to the best ones found
        Optionval=BestOptionval;
        Optdecision=Bestdecision;
        FdsumP=BestFdsumP;
        Fd0=BestFd0;
        V0=BestV0;
        Fdexpected=BestFdexpected;
        FixStoprule=BestFixStoprule;
    end
end
runtime=toc;
        

    