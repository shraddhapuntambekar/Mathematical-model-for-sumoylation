%% Mathematical model for sumoylation system - Parameter Variation
%%% Ref: Puntambekar, S. S., Nyayanit, D., Saxena, P., & Gadgil, C. J. (2016). Identification of Unintuitive Features of Sumoylation through Mathematical Modeling. Journal of Biological Chemistry, 291(18), 9458-9468.

%%% Studying the effect of parameter variation on the steady states of the
%%% simplest system (open system with no sumoylation of E2 and the second target)

%%% Parameters varied in the system:
%%% (1) SENP formation  (param.k(25,1))
%%% (2) binding of sumoE2 to the target  (param.k(14,1))
%%% (3) parameters for degradation of sumoE1 and sumoE20 (param.k(27,1), param.k(28,1))
clc; clear all; close all;

paravar = [1e-3 1e-2 1e-1  1 1e1  1e2 1e3];  %%% for k25, k14
interdeg = [0 0.001 0.01 0.1 1];            %%% for k27, k28

%%
getnames;
param=getparams_real;
t0=0; tf=1e2;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
options_fsol = optimset('Tolfun',eps,'MaxFunEvals',1e12,'MaxIter',1e12 );

%%% Simulating the simplest system
e2istar = 0; systemhas2tars = 0; isclosedsys=0;
if e2istar == 0, param.k(13:14,2)=0; param.k(16:19,:)=0; param.k(29,:)=0; end
if systemhas2tars == 0 , param.k(20:24,:)=0; end

y0=getinitial(param);

%% Parameter Variation
k14varied = param.k(14,1)*paravar;
k25varied = param.k(25,1)*paravar;
se2deg=param.k(27,1)*interdeg;
se1deg=param.k(28,1)*interdeg;

ck14 = size(paravar ,2);  ck2728=size(interdeg,2); csenp=size(k25varied,2);

t11_all=zeros(ck14,ck2728,csenp);
exitflagall=zeros(ck14,ck2728,csenp);

for ctrk14 =1:ck14
    param.k(14,1)=k14varied(ctrk14);
    
    for ctrk2728 = 1:ck2728
        param.k(27,1)=se2deg(ctrk2728);       %% para for sumoE20 degradation
        param.k(28,1)=se1deg(ctrk2728);       %% para for sumoE1 degradation
        
        for ctrsenp = 1:csenp
            param.k(25,1)=k25varied(ctrsenp);
            
            err =1;
            while err > 1e-8
                [t,y]=ode15s(@getrate,[t0 tf],y0,options,param);
                
                err = max(abs(y(end,:) - y0));
                y0=y(end,:);
                t0=tf;
                tf=tf+100;
            end
            
            %%% fsolve
            y0_forfsolve=y(end,:);
            for ctr1=1:param.n_species
                if (y0_forfsolve(ctr1)<0)
                    y0_forfsolve(ctr1)= 0;
                end
            end
            
            [x,fval,exitflag]= fsolve(@getrate_forfsolve,y0_forfsolve,options_fsol,param);
            
            t11_all(ck14,ck2728,csenp) = x(nt11);
            exitflagall(ck14,ck2728,csenp) = exitflag;
            
        end
    end
end

%%% save  k14_k25_k2829varied_open_e2nottar_1tar.mat


