%% May 2014

function param = getparams_All1

    param.n_species =18;  %%%(13 + 5 complexes)
    param.n_rxn =  42;%%30;  %%39;  %%% (27 + 12 )
    
    param.maxparaminrxn = 7; 
%     param.k(1:param.n_rxn,1:param.maxparaminrxn)=0; %matrix of 
param.k(1:31,1:param.maxparaminrxn)=0;
    
    param.cref= 1; %% 1 microMolar
    param.tref= 1;%60;%3e3;%60; %% 1 sec
     
    %%   
    param.k(1,1)=1*param.tref/param.cref;       %% preSumo formation
    param.k(2,1)=0.1*param.tref;                          %% preSumo degradation
    param.k(3,1)=0.1*param.tref;%%0.001*param.tref;                         %% degradation of SUMO
    param.k(4,1)=1*param.tref/param.cref;      %% E1 formation
    param.k(5,1)=0.1*param.tref;                         %% E1 degradation 
    param.k(6,1)=1*param.tref/param.cref;       %% E20 formation 
    param.k(7,1)=0.1*param.tref;                     %% E20 degradation 
    param.k(8,1)=0.525*param.tref/param.cref;  %% T10 formation  
    param.k(9,1)= 0.1*param.tref;                     %% T10 degradation 
    param.k(10,1)=0.1*param.tref;                          %% T11 degradation 
    
%     param.k(11,1)=1/param.cref;  %27.9/param.cref;   %% (27.9 +- 3.7 microM)  %% Km %% SUMO maturation - PROCESSING
%     param.k(11,2)= 1*param.tref/param.cref; %%3.6e-3*param.tref/param.cref;   %% Vmax (enzyme total of SENP2 is 5 nM)
%     param.k(11,3)=1*param.tref;%%0.72*param.tref;                       %%(0.72+-1.5 persec) %% Kcat
%     param.k(11,4)=param.k(11,3)/param.k(11,1); %%kcat/km (0.72/27.9 = 0.0258)
    param.k(11,5)=1;%0.2;%1*param.tref*param.cref;  %%% kf (Unit: perMole persecond)
    param.k(11,6)=0.1;%0.1*param.tref;  %%% kb (Unit:  persecond)
    param.k(11,7)=1;%0.2;%1*param.tref;  %%% kcat (Unit: persecond)
    
     %%%  (Kcat/Km = 0.78+-0.49/ 0.17+-0.096 per microMolar per sec) for  E1..rate constant from linear part of MM kinetics    
    param.k(12,1)=1*param.cref*param.tref;   %%4.5882%%kcat/km unit: 1/M.sec  %% sumo transfer to cys of E1 - ACTIVATION
    
     %%% (Kcat/Km = 1.06 +-0.05 per microMolar per sec)
    param.k(13,1)=1*param.tref*param.cref;     %1.06% sumo transfer from cys of E1 to cys of E20 - CONJUGATION to form product sumo~E20
    param.k(13,2)=1;%10*param.tref*param.cref;      %1.06% sumo transfer from cys of E1 to cys of E21 - CONJUGATION to form product sumo~E21
    
    %%% GST-RanGap1 without E3 (Kcat/Km = 0.66+-0.14/ 2.9+-0.95 = 0.2276 permicroMolar per sec)
    param.k(14,1)=1*param.tref*param.cref;     %0.2276% sumo transfer from cys of E20 to cys of T10 - LIGATION by E20 
    param.k(14,2)=1*param.tref*param.cref;   % 0.2276%VARIED%%% sumo transfer from cys of E21 to cys of T10 - LIGATION by E21
    
    %%% DECONJUGATION of T11 by SENP   
%     param.k(15,1)=1/param.cref;   %% Km = 33.1+- 3.2 MicroMolar %% Km 
%     param.k(15,2)=1*param.tref/param.cref;    %% Vmax %% Vmax = 25.2e-3 MicroMolar/sec
%     param.k(15,3)=0.5*param.tref;    %% Kcat %% Kcat =50.4 +-2.2 s-1
%     param.k(15,4)=param.k(15,3)/param.k(15,1); %%kcat/km ( 50.4/33.1 = 1.5227)
    param.k(15,5)=1*param.tref*param.cref;  %%% kf (Unit: perMole persecond)
    param.k(15,6)=0.1*param.tref;  %%% kb (Unit:  persecond)
    param.k(15,7)=1*param.tref;  %%% kcat (Unit: persecond)
    
    %% AUTOSUMOYLATION of E2 to form product E21 (E20*sumo)
    %%%ASSUMPTION (1) Will consider to have same parameter as that of conjugation reaction
    param.k(16,1)=1*param.tref*param.cref;                       %1.06% sumo transfer from cys of E1 to Lys of E20 - 
    % following parameters are if the E2 is modified as in Ubiquitination
%     param.k(16,2)=1;%0.1*param.tref*param.cref;                         %1.06% sumo transfer from cys of E20(sumo~E20) to lys of E20 - Ubiquitination like modification
%     param.k(16,3)=1;%0.1*param.tref*param.cref;                         %1.06% sumo transfer from cys of E21(sumo~E21) to lys of E20 - Ubiquitination like modification
    
    %%%ASSUMPTION (2) Deconjugation activity of SENP for E21 and sumoE21 is same as that for T11
     %%% DECONJUGATION of E21 by SENP  - considering it to be same as that for T11
%     param.k(17,1)=1/param.cref;   %% Km = 33.1+- 3.2 MicroMolar %% Km %% DECONJUGATION of T11 by SENP   
%     param.k(17,2)=1*param.tref/param.cref;    %% Vmax %% Vmax = 25.2e-3 MicroMolar/sec
%     param.k(17,3)=0.5*param.tref;    %% Kcat %% Kcat =50.4 +-2.2 s-1
%     param.k(17,4)=param.k(17,3)/param.k(17,1); %%kcat/km
    param.k(17,5)=0.2;%0;%0.2;%*param.tref*param.cref;  %%% kf (Unit: perMole persecond)
    param.k(17,6)=0.05;%0;%0.05;%0.01*param.tref;  %%% kb (Unit:  persecond)
    param.k(17,7)=0.2;%0;%0.2;%0.5*param.tref;  %%% kcat (Unit: persecond)

    %%% DECONJUGATION of sumoE21 by SENP  - considering it to be same as that for T11
%     param.k(18,1)=1/param.cref;   %% Km = 33.1+- 3.2 MicroMolar %% Km %% DECONJUGATION of T11 by SENP   
%     param.k(18,2)=1*param.tref/param.cref;    %% Vmax %% Vmax = 25.2e-3 MicroMolar/sec
%     param.k(18,3)=0.5*param.tref;    %% Kcat %% Kcat =50.4 +-2.2 s-1
%     param.k(18,4)=param.k(18,3)/param.k(18,1); %%kcat/km
    param.k(18,5)=0.2;%0;%0.2;%1*param.tref*param.cref;  %%% kf (Unit: perMole persecond)
    param.k(18,6)=0.05;%0;%0.05;%0.01*param.tref;  %%% kb (Unit:  persecond)
    param.k(18,7)=0.2;%0;%0.2;%0.5*param.tref;  %%% kcat (Unit: persecond)
    
    param.k(19,1)=0.1*param.tref;                  %% degradation of E21
    
    param.k(20,1)=0; %% k(14,2) - in DImpals code - UNUSED TO KEEP UNIFORMITY IN PARAMETER NUMBERING
    param.k(21,1)=0;%% k(13,2) - in DImpals code - UNUSED TO KEEP UNIFORMITY IN PARAMETER NUMBERING
    
    %% SECOND TARGET
    param.k(22,1)=1*param.tref/param.cref;  %% T20 formation  
    param.k(23,1)= 0.1*param.tref;                     %% T20 degradation 
    param.k(24,1)=0.1*param.tref;                          %% T21 degradation 
    
    %%% GST-RanGap1 without E3 (Kcat/Km = 0.66+-0.14/ 2.9+-0.95 = 0.2276 permicroMolar per sec)
    param.k(25,1)=1*param.tref*param.cref;     %% sumo transfer from cys of E20 to cys of T10 - LIGATION by E20 
    param.k(25,2)=1*param.tref*param.cref;    %VARIED%%% sumo transfer from cys of E21 to cys of T10 - LIGATION by E21
    
    %%% DECONJUGATION of T21 by SENP   
%     param.k(26,1)=1/param.cref;   %% Km = 33.1+- 3.2 MicroMolar %% Km 
%     param.k(26,2)=1*param.tref/param.cref;    %% Vmax %% Vmax = 25.2e-3 MicroMolar/sec
%     param.k(26,3)=0.25*param.tref;    %% Kcat %% Kcat =50.4 +-2.2 s-1
%     param.k(26,4)=param.k(26,3)/param.k(26,1); %%kcat/km ( 50.4/33.1 = 1.5227)
    param.k(26,5)=0.1*param.tref*param.cref;  %%% kf (Unit: perMole persecond)
    param.k(26,6)=0.025*param.tref;  %%% kb (Unit:  persecond)
    param.k(26,7)=0.1*param.tref;  %%% kcat (Unit: persecond)
    
    %%
    param.k(27,1)=1*param.tref/param.cref;      %% SENP formation
    param.k(28,1)=0.1*param.tref;                      %% SENP degradation
    
    param.k(29,1)=0.01*param.tref;%0;%0.1*param.tref;        %% sumoE20 degradation
    param.k(30,1)=0.01*param.tref;%0;%0.1*param.tref;        %% sumoE1 degradation
  param.k(31,1)=0.01*param.tref;%0;%0.1*param.tref;        %% sumoE21 degradation
return;
