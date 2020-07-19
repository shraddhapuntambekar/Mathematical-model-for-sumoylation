function dy = getrate(t,y,param)
getnames;
s = zeros(param.n_species,param.n_rxn);
r = zeros(param.n_rxn,1);

yc = num2cell(y);
[presumo,sumo,senp,e10,se10,e20,e21,se20,se21,t10,t11,t20,t21,presumosenp,t11senp,t21senp,e21senp,sumoe21senp]=deal(yc{:});

%%%Pool->presumo
r(1,1) = param.k(1,1);
s(npresumo,1) = 1; 

%%%presumo->Pool
r(2,1) = param.k(2,1)*presumo;
s(npresumo,2) = -1; 

%%% sumo->Pool
r(3,1) = param.k(3,1)*sumo;
s(nsumo,3) = -1; 

%%% Pool->e10
r(4,1) = param.k(4,1);
s(ne10,4) = 1; 

%%% e10->Pool
r(5,1) = param.k(5,1)*e10;
s(ne10,5) = -1; 

%%% Pool->e20
r(6,1) = param.k(6,1);
s(ne20,6) = 1; 

%%% e20->Pool
r(7,1) = param.k(7,1)*e20;
s(ne20,7) = -1; 

%%% Pool->tar
r(8,1) = param.k(8,1);
s(nt10,8) = 1; 

%%% t10->Pool
r(9,1) = param.k(9,1)*t10;
s(nt10,9) = -1; 

%%% t11->Pool
r(10,1) = param.k(10,1)*t11;
s(nt11,10) = -1;

%%
%%% presumo + senp-> [presumosenp]  --Preprocessing
r(11,1) = param.k(11,1)*presumo*senp;
s(npresumo,11) = -1 ;s(nsenp,11)= -1; s(npresumosenp,11)= 1; 

%%% sumo +e10 ->se10 -- Activation
r(12,1) = param.k(12,1)*sumo*e10;
s(nsumo,12) = -1 ; s(ne10,12)= -1; s(nse10,12)= 1;  

%%% se10+e20->se20+e10 - Conjugation
r(13,1) = param.k(13,1)*se10*e20;
s(nse10,13) = -1; s(ne20,13)= -1; s(nse20,13)= 1; s(ne10,13)=1; 

%se20+t10->t11+e20 - Ligation of tar1
r(14,1) = param.k(14,1)*se20*t10;
s(nse20,14) = -1; s(nt10,14)= -1; s(nt11,14)=1; s(ne20,14)=1; 

%%% t11+senp->[t11senp] - Deconjugation of T11 by SENP
r(15,1) = param.k(15,1)*t11*senp;
s(nt11,15) = -1; s(nsenp,15)=-1; s(nt11senp,15)= 1; 

%% when e2 is target
%se10+e20->e21+e10 - Modification of E2 by sumo~E1
r(16,1) = param.k(16,1)*se10*e20;
s(nse10,16) = -1; s(ne20,16) = -1; s(ne21,16) = 1; s(ne10,16) = 1; 
%%% Modification of E2 like in Ubiquitination se20+e20->e21+e20 
% % % r(16,1) = param.k(16,2)*se10*e20;
% % % s(nse10,16) = -1; s(ne20,16) = -1; s(ne21,16) = 1; s(ne10,16) = 1; 

%%% e21+senp-> [e21senp] -- Deconjugation of e21 by SENP
r(17,1) = param.k(17,1)*e21*senp;
s(ne21,17) = -1; s(nsenp,17)= -1; s(ne21senp,17)= 1; 

%se21+senp-> [se21senp] -- Deconjugation of se21 by SENP
r(18,1) = param.k(18,1)*se21*senp;
s(nse21,18) = -1; s(nsenp,18)= -1; s(nsumoe21senp,18)= 1; 

%e21->Pool
r(19,1) = param.k(19,1)*e21;
s(ne21,19) = -1; 

%%
%se21+t10->t11+e21 - Ligation of target by modified E2
r(20,1) = param.k(14,2)*se21*t10;
s(nse21,20) = -1; s(nt10,20)= -1; s(nt11,20)=1; s(ne21,20)=1;

%se10+e21->se21+e10 - Conjugation of modified E2
r(21,1) = param.k(16,1)*se10*e21;
s(nse10,21) = -1; s(ne21,21) = -1; s(nse21,21) = 1; s(ne10,21) = 1; 

%% second target in the system
%%% Pool->tar
r(22,1) = param.k(20,1);
s(nt20,22) = 1; 

%%% t20->Pool
r(23,1) = param.k(21,1)*t20;
s(nt20,23) = -1; 

%%% t21->Pool
r(24,1) = param.k(22,1)*t21;
s(nt21,24) = -1; 

%%
%%% se20+t20->e20+t21 
r(25,1) = param.k(23,1)*se20*t20;
s(nse20,25) = -1; s(nt20,25)= -1; s(ne20,25)= 1; s(nt21,25)=1; 
 
%%% se21+t20->e21+t21 
r(26,1) = param.k(23,2)*se21*t20;
s(nse21,26) = -1; s(nt20,26)=-1; s(ne21,26)= 1; s(nt21,26)=1; 

%%% t21+senp-> [t21senp]
r(27,1) = param.k(24,1)*t21*senp;
s(nt21,27) = -1; s(nsenp,27)=-1; s(nt21senp,27)= 1; 

%% added reactions of complexes with SENP
%%%  [presumosenp]--> presumo + senp  --Preprocessing
r(28,1) = param.k(11,2)*presumosenp;
s(npresumo,28) = 1 ;s(nsenp,28)= 1; s(npresumosenp,28)= -1; 

%%%  [presumosenp]--> sumo + senp  --Preprocessing
r(29,1) = param.k(11,3)*presumosenp;
s(nsumo,29) = 1 ;s(nsenp,29)= 1; s(npresumosenp,29)= -1; 

%%% [t11senp] -> t11+senp - Deconjugation of T11 by SENP
r(30,1) = param.k(15,2)*t11senp;
s(nt11,30) = 1; s(nsenp,30)=1; s(nt11senp,30)= -1; 

%%% [t11senp] -> t10+senp + sumo - Deconjugation of T11 by SENP
r(31,1) = param.k(15,3)*t11senp;
s(nt10,31) = 1; s(nsenp,31)=1; s(nsumo,31)=1; s(nt11senp,31)= -1;

%%% [e21senp] --> e21+senp -- Deconjugation of e21 by SENP
r(32,1) = param.k(17,2)*e21senp;
s(ne21,32) = 1; s(nsenp,32)= 1; s(ne21senp,32)= -1; 

%%% [e21senp] --> e20+senp+sumo -- Deconjugation of e21 by SENP
r(33,1) = param.k(17,3)*e21senp;
s(ne20,33) = 1; s(nsenp,33)= 1; s(ne21senp,33)= -1; s(nsumo,33)=1;

%se21+senp-> [se21senp] -- Deconjugation of se21 by SENP
r(34,1) = param.k(18,2)*sumoe21senp;
s(nse21,34) = 1; s(nsenp,34)= 1; s(nsumoe21senp,34)= -1; 

%%% [se21senp] --> se20+senp + sumo -- Deconjugation of se21 by SENP
r(35,1) = param.k(18,3)*sumoe21senp;
s(nse20,35) = 1; s(nsenp,35)= 1; s(nsumoe21senp,35)= -1; s(nsumo,35)=1;

%%%  [t21senp] --> t21+senp
r(36,1) = param.k(24,2)*t21senp;
s(nt21,36) = 1; s(nsenp,36)=1; s(nt21senp,36)= -1; 

%%%  [t21senp] --> t20+senp + sumo
r(37,1) = param.k(24,3)*t21senp;
s(nt20,37) = 1; s(nsenp,37)=1; s(nt21senp,37)= -1; s(nsumo,37)=1;

%%
%%%Pool->senp
r(38,1) = param.k(25,1);
s(nsenp,38) = 1; 

%%%senp->Pool
r(39,1) = param.k(26,1)*senp;
s(nsenp,39) = -1; 

%%% sumoE20 ->Pool
r(40,1) = param.k(27,1)*se20;
s(nse20,40) = -1;

%%% sumoE1 ->Pool
r(41,1) = param.k(28,1)*se10;
s(nse10,41) = -1;

%%% sumoE21 ->Pool
r(42,1) = param.k(29,1)*se21;
s(nse21,42) = -1;

dy =s*r; 
                                                        

 