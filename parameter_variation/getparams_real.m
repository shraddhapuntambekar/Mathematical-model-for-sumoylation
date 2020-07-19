function param = getparams_real

    param.n_species =18;  %%%(13 + 5 complexes)
    param.n_rxn =  42;
    
    param.maxparaminrxn = 3; 
    param.k(1:29,1:param.maxparaminrxn)=0;
 
    param.cref= 0.01;    %%  microMolar
    param.tref= 1;       %%  hr

%%   
    param.k(1,1)=6.286e-03*param.tref/param.cref;      %%% preSumo formation (microM/hr)
    param.k(2,1)=(log(2)/(13.57))*param.tref;          %%% preSumo degradation t1/2=13.57 h (assumed)
    param.k(3,1)=(log(2)/(13.57))*param.tref;          %%% degradation of SUMO t1/2=13.57 h
    
    param.k(4,1)=2.426e-2*param.tref/param.cref;       %%% E1a formation (microM/hr)
    param.k(5,1)=(log(2)/(69.79))*param.tref;          %%% E1a degradation  t1/2=69.79 h
    
% param.k(4,1)=6.3010322e-3*param.tref/param.cref;     %%% E1b formation (microM/hr)
% param.k(5,1)=(log(2)/(102.92))*param.tref;           %%% E1b degradation 
      
    param.k(6,1)=4.416e-2*param.tref/param.cref;       %%% E20 formation (microM/hr)
    param.k(7,1)=(log(2)/(72.85))*param.tref;          %%% E20 degradation 
    
    param.k(8,1)=4.591e-2*param.tref/param.cref;       %%% T10 formation  (microM/hr)
    param.k(9,1)= (log(2)/(15.53))*param.tref;         %%% T10 degradation 
    param.k(10,1)=(log(2)/(15.53))*param.tref;         %%% T11 degradation (Assuming)
    
    %% 
    
    param.k(11,1)=1*3600*param.tref*param.cref;        %%% kf (Unit: perMole persecond) (assuming kf=1) - %%% Preprocessing of presumo
    param.k(11,2)=27.18*3600*param.tref;               %%% kb (Unit:  persecond)  (kr= km*kf  - kcat) .. kr = 27.9 - 0.72  (km is in micromolar)
    param.k(11,3)=0.72*3600*param.tref;                %%% kcat (Unit: persecond)
    
     %%%  (Kcat/Km = 0.78+-0.49/ 0.17+-0.096 per microMolar per sec) for  E1..rate constant from linear part of MM kinetics    
    param.k(12,1)=4.5882*3600*param.cref*param.tref;   %%% kcat/km unit: 1/microM.sec  %% sumo transfer to cys of E1 - ACTIVATION
    
     %%% (Kcat/Km = 1.06 +-0.05 per microMolar per sec)
    param.k(13,1)=1.06*3600*param.tref*param.cref;     %%% sumo transfer from cys of E1 to cys of E20 - CONJUGATION to form product sumo~E20
    param.k(13,2)=1.06*3600*param.tref*param.cref;     %%% sumo transfer from cys of E1 to cys of E21 - CONJUGATION to form product sumo~E21
    
    %%% GST-RanGap1 without E3 (Kcat/Km = 0.66+-0.14/ 2.9+-0.95 = 0.2276 permicroMolar per sec)
    param.k(14,1)=0.2276*3600*param.tref*param.cref;   %%% sumo transfer from cys of E20 to cys of T10 - LIGATION by E20 
    param.k(14,2)=0.2276*3600*param.tref*param.cref;   %%% sumo transfer from cys of E21 to cys of T10 - LIGATION by E21
    
    %%% DECONJUGATION of T11 by SENP   
    param.k(15,1)=10*3600*param.tref*param.cref;       %%% kf (Unit: perMole persecond)
    param.k(15,2)=280.6*3600*param.tref;               %%% kb (Unit:  persecond)
    param.k(15,3)=50.4*3600*param.tref;                %%% kcat (Unit: persecond)
    
    %% AUTOSUMOYLATION of E2 to form product E21 (E20*sumo)
    %%%ASSUMPTION (1) Will consider to have same parameter as that of conjugation reaction
    param.k(16,1)=1.06*3600*param.tref*param.cref;     %%% sumo transfer from cys of E1 to Lys of E20 - 
%     % following parameters are if the E2 is modified as in Ubiquitination
%     param.k(16,2)=1.06*3600*param.tref*param.cref;   %%% sumo transfer from cys of E20(sumo~E20) to lys of E20 - Ubiquitination like modification
%     param.k(16,3)=1.06*3600*param.tref*param.cref;   %%% sumo transfer from cys of E21(sumo~E21) to lys of E20 - Ubiquitination like modification
    
    %%%ASSUMPTION (2) Deconjugation activity of SENP for E21 and sumoE21 is same as that for T11
     %%% DECONJUGATION of E21 by SENP  - considering it to be same as that for T11
    param.k(17,1)=10*3600*param.tref*param.cref;       %%% kf (Unit: perMole persecond)
    param.k(17,2)=280.6*3600*param.tref;               %%% kb (Unit:  persecond)
    param.k(17,3)=50.4*3600*param.tref;                %%% kcat (Unit: persecond)

    %%% DECONJUGATION of sumoE21 by SENP  - considering it to be same as that for T11
    param.k(18,1)=10*3600*param.tref*param.cref;       %%% kf (Unit: perMole persecond)
    param.k(18,2)=280.6*3600*param.tref;               %%% kb (Unit:  persecond)
    param.k(18,3)=50.4*3600*param.tref;                %%% kcat (Unit: persecond)
    
    param.k(19,1)=(log(2)/(72.85))*param.tref;         %%% degradation of E21
       
    %% SECOND TARGET
    param.k(20,1)=4.591e-2*param.tref/param.cref;             %%% T20 formation  
    param.k(21,1)=(log(2)/(15.53))*param.tref;                     %%% T20 degradation 
    param.k(22,1)=(log(2)/(15.53))*param.tref;                      %%% T21 degradation 
    
    %%% GST-RanGap1 without E3 (Kcat/Km = 0.66+-0.14/ 2.9+-0.95 = 0.2276 permicroMolar per sec)
    param.k(23,1)=0.2276*3600*param.tref*param.cref;             %%% sumo transfer from cys of E20 to cys of T10 - LIGATION by E20 
    param.k(23,2)=0.2276*3600*param.tref*param.cref;             %%% sumo transfer from cys of E21 to cys of T10 - LIGATION by E21
    
    %%% DECONJUGATION of T21 by SENP   
    param.k(24,1)=10*3600*param.tref*param.cref;           %%% kf (Unit: perMole persecond)
    param.k(24,2)=280.6*3600*param.tref;                    %%% kb (Unit:  persecond)
    param.k(24,3)=50.4*3600*param.tref;                     %%% kcat (Unit: persecond)
    
    %%
    param.k(25,1)=3.212e-4*param.tref/param.cref;      %%% SENP formation (microM/hr)
    param.k(26,1)=(log(2)/(20.07))*param.tref;         %%% SENP degradation
       
    param.k(27,1)=(log(2)/(72.85))*param.tref;         %%% sumoE20 degradation
    param.k(28,1)=(log(2)/(69.79))*param.tref;         %%% sumoE1 degradation
    param.k(29,1)=(log(2)/(72.85))*param.tref;         %%% sumoE21 degradation
return;
