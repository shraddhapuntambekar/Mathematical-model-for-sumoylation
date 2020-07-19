function param = getparams(numofSUMOontarget,numOFconstspecies)

    numOFconstrexn = 7;  %%%{(1n2)presumo form n degrade (3n4) sumo form n degrade (5)processing of presumo}
    %numOFvarrexn = 7;   %%%{(1)activation (2)conjugation (3)ligation of tars n e1 (4)ligation of e2 (5)deconjugation (6)tar formation (7)tar degradatn}
    cref=1;              %%%  microMolar
    tref=1;              %%%  sec
    %%
    param.numofSUMOontarget=numofSUMOontarget;
    param.Totrxn= numOFconstrexn + (param.numofSUMOontarget(2)+1) + ((param.numofSUMOontarget(2)+1) * (param.numofSUMOontarget(1)+1)) + ((param.numofSUMOontarget(1)+1) * sum(param.numofSUMOontarget(2:end))) + ((param.numofSUMOontarget(2)+1) * sum(param.numofSUMOontarget(1))) + sum(param.numofSUMOontarget) + max(size(param.numofSUMOontarget)) + (sum(param.numofSUMOontarget) + max(size(param.numofSUMOontarget))) ;
    Totsps= numOFconstspecies + sum(2*(param.numofSUMOontarget(1:2)+1)) + sum(param.numofSUMOontarget(3:end)+1);

    param.n_species = Totsps;
    param.numOFconstrexn=numOFconstrexn;
    param.numOFconstspecies=numOFconstspecies;
    param.maxparaminrxn = 3; % for other constants such as km
    param.k(1:param.numOFconstrexn,1:param.maxparaminrxn)=0;

    param.cref= cref;
    param.tref= tref;

    param.k(1,1)=1*param.tref/param.cref;       %% preSumo formation
    param.k(2,1)=1*param.tref;                  %% preSumo degradation

    param.k(3,1)=1*param.tref/param.cref;       %%  SUMO  formation
    param.k(4,1)=1*param.tref;                  %% degradation of SUMO
    
    param.k(5,1)=1*param.tref/param.cref;       %% E3 formation
    param.k(6,1)=1*param.tref;                  %% E3 degradation
     
%     param.k(7,1)= 1/param.cref;           %% Km %% SUMO maturation -Processing
%     param.k(7,2)=1*param.tref/param.cref; %% Vmax
    
    param.k(7,1)= 27.9/param.cref;   %% (27.9 +- 3.7 microM)  %% Km %% SUMO maturation -Processing
    param.k(7,2)=3.6e-3*param.tref/param.cref;    %% (3.6e-3 MicroM/sec)   %% Vmax
   
    param.ke3 = 1;   % rate for E3 activity

    %%
    param.k8_activation  = ones(1,param.numofSUMOontarget(2)+1);  %% k8_activation=[E10 E11 E12 . . E1n] Unit: k8 = per Mol per sec
    param.k8_activation(1,1) = 0.1560*param.tref*param.cref;  %%% Conidering maximal rate(Vmax) for E1   
    
    %                      sumoE10 sumoE11 sumoE12 . . sumoE1n
    % % %%k9_conjugation = [                                        E20
    % %                                                             E21
    % %                                                               .
    % %                                                           ] E2n
    param.k9_conjugation = ones((param.numofSUMOontarget(1)+1),(param.numofSUMOontarget(2)+1))*param.tref*param.cref; %%% units: per mole per sec
    param.k9_conjugation(1,:)=param.k9_conjugation(1,:)*0.0323e-6; %%% for E20
    param.k9_conjugation(2:end,:)=param.k9_conjugation(2:end,:); %%% sumoylates E2 can be made non enymatic rendering no effecct by setting parameter to 0

    %                      sumoE20 sumoE21 sumoE22 . . sumoE2n
    % %k10tare1_ligation = [                                        E10
    % %                                                            E11
    % %                                                               .
    % %                                                            E1(n-1)
    % %                                                            T10
    % %                                                          ] T1(n-1)
%     param.k10tare1_ligation = ones(sum(param.numofSUMOontarget(2:end)),(param.numofSUMOontarget(1)+1));
%     param.k10tare1_ligation(end-param.numofSUMOontarget(end)+1:end,:)=param.k10tare1_ligation(end-param.numofSUMOontarget(end)+1:end,:)*0; %%last target has slow reaction ligation activity
    
    if param.numofSUMOontarget(2) == 0
        param.k10e1_ligation = zeros(1,(param.numofSUMOontarget(1)+1))*param.tref*param.cref;
    else
        param.k10e1_ligation = ones(param.numofSUMOontarget(2),(param.numofSUMOontarget(1)+1))*param.tref*param.cref;
    end


    param.k10tar_ligation = ones(sum(param.numofSUMOontarget(3:end)),(param.numofSUMOontarget(1)+1))*param.tref*param.cref;
    %param.k10tar_ligation(end-param.numofSUMOontarget(end)+1:end,:)=param.k10tar_ligation(end-param.numofSUMOontarget(end)+1:end,:)*0.01; %%last target has slow reaction ligation activity
    param.k10tar_ligation(:,2)=param.k10tar_ligation(:,2); %can make sumoE21, sumoe22 inactive if parameter is set to 0

    param.k10tar_ligation(end-param.numofSUMOontarget(end)+1:end,:)=param.k10tar_ligation(end-param.numofSUMOontarget(end)+1:end,:)*180e-6; %% p53 180 M s-1
    
    
    
    %                       sumoE10 sumoE11 sumoE12 . . sumoE1n
    % %k11e2_ligation =    [                                        E20
    % %                                                            E21
    % %                                                               .
    % %                                                          ] E2(n-1)
    param.k11e2_ligation = ones(param.numofSUMOontarget(1),(param.numofSUMOontarget(2)+1))*param.tref*param.cref*0.1; 
    
% % %                            E10. .E1(n-1) E20. .E2(n-1) T10. .T1(n-1) 
% % %     %%% k12_deconjugation=[                                         ; %km
% % %                                                                     ] %Vmax 

% %     param.k12_deconjugation = ones(2,sum(param.numofSUMOontarget));  
% %     param.k12_deconjugation(2,:)=param.k12_deconjugation(2,:)*0; % Vmax

    param.k12_deconjugation = ones(2,sum(param.numofSUMOontarget));
    param.k12_deconjugation(1,:)=param.k12_deconjugation(1,:)/param.cref; %%% Km   %% non parameterize
    param.k12_deconjugation(2,:)=param.k12_deconjugation(2,:)*param.tref/param.cref; %%% Vmax %% non parameterize
    
     %%% For target1>> RanGap1
    param.k12_deconjugation(1,sum(param.numofSUMOontarget(1,1:2)) + 1)=33.1 ; %% Km = 33.1 MicroMolar
    param.k12_deconjugation(2,sum(param.numofSUMOontarget(1,1:2)) + 1)=25.2e-3; %% Vmax = 25.2e-3 MicroMolar/sec
    
    param.k12_deconjugation(2,:)=param.k12_deconjugation(2,:); % Vmax  (to make closed sys)
     
    param.k13_e2e1tarform = ones(1,size(param.numofSUMOontarget,2)); %%param.k13_e2e1tarform =[E20 E10 T10 T20..]

% % %                         E2  E1  T1  T2 
% % % %%% k14_e2e1tardegrade=[                 ; Ts  
% % %                                          ] Tss    
    param.k14_e2e1tardegrade = ones(max(param.numofSUMOontarget) + 1,size(param.numofSUMOontarget,2));
return;
