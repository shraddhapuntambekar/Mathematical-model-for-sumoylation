function dy = getrate(t,y,param)

getname;
s = zeros(param.n_species,param.Totrxn);
r = zeros(param.Totrxn,1);

yc = num2cell(y);
[presumo,sumo,senp,e3]=deal(yc{1:param.numOFconstspecies});

 %%
        %%%Reaction 1::  Pool->presumo  %Presumo formation 
        r(1,1) = param.k(1,1);
        s(npresumo,1) = 1;

        %%%Reaction 2::  presumo -> pool %%Presumo degradation 
        r(2,1) = param.k(2,1)*presumo;
        s(npresumo,2) = -1;

        %%%Reaction 3:: pool -> sumo  %%Sumo formation
        r(3,1) = param.k(3,1);
        s(nsumo,3) = 1;
        
        %%%Reaction 4:: sumo ->pool    %% Sumo degradation
        r(4,1) = param.k(4,1)*sumo;
        s(nsumo,4) = -1;
         
        %%%Reaction 5::  Pool->E3   % E3 formation 
        r(5,1) = param.k(5,1);
        s(ne3,5) = 1;

        %%%Reaction 6:: E3 ->pool   %% E3 degradation 
        r(6,1) = param.k(6,1)*e3;
        s(ne3,6) = -1;      
        
        %%%Reaction 7:: presumo + senp-> sumo   %% Processing
        r(7,1) = (param.k(7,2)*presumo)/(param.k(7,1)+presumo);
        s(npresumo,7) = -1 ;s(nsumo,7)= 1;

  %%
        %%%Reaction 8:: sumo + e1 ->sumoe1(se1)  %% Activation
        indrex8 = param.numOFconstrexn;
        indE1 = param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1));
        indsumoE1 =  param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1)) + param.numofSUMOontarget(2) + 1;
        for ctr8 = 1:param.numofSUMOontarget(2)+1
            
            r(indrex8 + ctr8,1) = param.k8_activation(1,ctr8)* sumo * y(indE1 + ctr8);       % r(8,1) = param.k(8,1)*sumo*e1;
            s(nsumo,indrex8 + ctr8) = -1 ;           % s(nsumo,8) = -1 ; 
            s(indE1 + ctr8,indrex8 + ctr8)= -1;      % s(ne1,8)= -1;
            s(indsumoE1 + ctr8 , indrex8 + ctr8)= 1; % s(nse1,8)= 1;
            
        end
        
         %%%Reaction 9:: se1+e2->se2+e1    %% Conjugation
         indrex9 = indrex8 + ctr8;
         indE2 = param.numOFconstspecies;
         indsumoE2 = param.numOFconstspecies + param.numofSUMOontarget(1) + 1;
         tmp=0;
         for ctr9_e2 = 1: param.numofSUMOontarget(1) + 1
             for ctr9_se1 =1:param.numofSUMOontarget(2)+1
                  tmp=tmp+1;
                  
                  r(indrex9 + tmp,1) = (param.k9_conjugation(ctr9_e2,ctr9_se1)*y(indE2 + ctr9_e2)*y(indsumoE1 + ctr9_se1)) + (param.ke3 * e3*y(indE2 + ctr9_e2)*y(indsumoE1 + ctr9_se1));  % r(9,1) = param.k(9,1)*e2*se1 + (ke3*e3);
                  s(indsumoE1 + ctr9_se1 , indrex9 + tmp) = -1;  % s(nse1,9)=-1; 
                  s(indE2 + ctr9_e2 ,indrex9 + tmp)= -1;   % s(ne2,9)= -1; 
                  s(indE1 + ctr9_se1 , indrex9 + tmp) = 1;   % s(ne1,9) = 1; 
                  s(indsumoE2 + ctr9_e2 , indrex9 + tmp)=1;  % s(nse2,9)=1;
                  
             end
         end
  
         %%%Reaction 10:: se2 + e10/ target -> e11/target*sumo + e2    %% Ligation of E1 and Targets
         tmpr=1;
         tmp10=0;
         indrex10 = indrex9 + tmp;
         indtar =  param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1)) + (2*(param.numofSUMOontarget(2) + 1));
         for ctr10_e1tar = 2: size(param.numofSUMOontarget,2)
            for ctr10_sumoONe1tar = 1: param.numofSUMOontarget(ctr10_e1tar)
                 for ctr10_se2 =1:param.numofSUMOontarget(1)+1
                     tmp10=tmp10+1;
                                          
                     if ctr10_e1tar == 2
                           r(indrex10 + tmp10,1) = (param.k10e1_ligation(ctr10_sumoONe1tar,ctr10_se2) * y(indsumoE2 + ctr10_se2) * y(indE1 + ctr10_sumoONe1tar) ) + (param.ke3 * e3* y(indsumoE2 + ctr10_se2) * y(indE1 + ctr10_sumoONe1tar)); % r(10,1) = param.k(10,1)*se2*tar + (ke3*e3);
                           s(indE1 + ctr10_sumoONe1tar,indrex10 + tmp10)= -1; % s(ne10/tar,10)= -1;
                           s(indE1 + ctr10_sumoONe1tar + 1, indrex10 + tmp10 ) = 1; % s(e11/target*sumo,10) = 1;
                         
                     else
                           r(indrex10 + tmp10,1) = (param.k10tar_ligation(tmpr,ctr10_se2) * y(indsumoE2 + ctr10_se2) * y(indtar + ctr10_sumoONe1tar) ) + (param.ke3 * e3* y(indsumoE2 + ctr10_se2) * y(indtar + ctr10_sumoONe1tar)); % r(10,1) = param.k(10,1)*se2*tar + (ke3*e3);
                           s( indtar + ctr10_sumoONe1tar,indrex10 + tmp10)= -1; % s(ne10/tar,10)= -1;
                           s( indtar + ctr10_sumoONe1tar + 1, indrex10 + tmp10 ) = 1; % s(e11/target*sumo,10) = 1;
                           
                     end
                     

                     
                     s(indsumoE2 + ctr10_se2 , indrex10 + tmp10) = -1;  % s(nse2,10)=-1;
                     s(indE2 + ctr10_se2 ,indrex10 + tmp10) = 1;    % s(e2,10)=1;
                         
                 end
            end
             
            if ctr10_e1tar > 2
                indtar = indtar + param.numofSUMOontarget(ctr10_e1tar) + 1;
                tmpr=tmpr+1;
            end

         end
         
         %%%Reaction 11:: se10 + e20 -> e21 + e10  %%%% Ligation of E2 - Transfer of sumo on lys of e2 from cys of e1
         tmp11 = 0;
         indrex11 = indrex10 + tmp10;
         for ctr11_e2 = 1: param.numofSUMOontarget(1)
             for ctr11_se1 =1:param.numofSUMOontarget(2)+1
                 tmp11 = tmp11 + 1;
 
                   r(indrex11 + tmp11,1) =  (param.k11e2_ligation(ctr11_e2 , ctr11_se1) * y(indsumoE1 + ctr11_se1) * y(indE2 + ctr11_e2)) + (param.ke3 * e3* y(indsumoE1 + ctr11_se1) * y(indE2 + ctr11_e2));  % r(11,1) = param.k(11,1)*se1*e2 + (ke3*e3);
                   s(indsumoE1 + ctr11_se1,indrex11 + tmp11)= -1;  % s(nse1,11)= -1;
                   s(indE1 + ctr11_se1,indrex11 + tmp11) = 1;  % s(e1,11) = 1;
                   s(indE2 + ctr11_e2 ,indrex11 + tmp11) =-1;  % s(ne20,11)=-1;
                   s(indE2 + ctr11_e2 + 1,indrex11 + tmp11)=1;  % s(e21,11)=1;
             end
         end
         
        %%%Reaction 12:: Tar*sumo -(senp)-> Tar + sumo  %%%% Deconjugation
        tmp12 = 0;
        indrex12 = indrex11 + tmp11;
        indtar =  param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1)) + (2*(param.numofSUMOontarget(2) + 1));
        for ctr12_ovrtar = 1:size(param.numofSUMOontarget,2)
            for ctr12_sumoONtar = 1:param.numofSUMOontarget(ctr12_ovrtar)
                tmp12 = tmp12 + 1;
                
                if ctr12_ovrtar == 1   %% E2
                    r(indrex12 + tmp12,1) = (param.k12_deconjugation(2,tmp12)*y(indE2 + ctr12_sumoONtar + 1))/(param.k12_deconjugation(1,tmp12) + y(indE2 + ctr12_sumoONtar + 1)); % r(12,1) = (param.k(12,2)*tar-sumo)/(param.k(12,1)+tar-sumo);
                    s(indE2 + ctr12_sumoONtar + 1 ,indrex12 + tmp12) = -1 ; % s(Tar*sumo,12) = -1 ;
                    s(indE2 + ctr12_sumoONtar,indrex12 + tmp12)= 1; % s(Tar,12)= 1;
                elseif ctr12_ovrtar == 2  %% E1
                    r(indrex12 + tmp12,1) = (param.k12_deconjugation(2,tmp12)*y(indE1 + ctr12_sumoONtar + 1))/(param.k12_deconjugation(1,tmp12) + y(indE1 + ctr12_sumoONtar + 1)); % r(12,1) = (param.k(12,2)*tar-sumo)/(param.k(12,1)+tar-sumo);
                    s(indE1 + ctr12_sumoONtar + 1 ,indrex12 + tmp12) = -1 ; % s(Tar*sumo,12) = -1 ;
                    s(indE1 + ctr12_sumoONtar ,indrex12 + tmp12) = 1 ; % s(Tar,12)= 1;
                else  %% Targets
                    r(indrex12 + tmp12,1) = (param.k12_deconjugation(2,tmp12)*y(indtar + ctr12_sumoONtar + 1))/(param.k12_deconjugation(1,tmp12) + y(indtar + ctr12_sumoONtar + 1)); % r(12,1) = (param.k(12,2)*tar-sumo)/(param.k(12,1)+tar-sumo);
                    s(indtar + ctr12_sumoONtar + 1 ,indrex12 + tmp12) = -1 ; % s(Tar*sumo,12) = -1 ;
                    s(indtar + ctr12_sumoONtar ,indrex12 + tmp12) = 1 ;% s(Tar,12)= 1;
                end
                
                s(nsumo,indrex12 + tmp12) = 1 ; % s(sumo,12) = 1 ;
            end
            
                if ctr12_ovrtar > 2
                    indtar = indtar + param.numofSUMOontarget(ctr12_ovrtar) + 1;
                end
            
        end
        
       %%%% Target formation and Target Degradation:: 
       indrex13 = indrex12 + tmp12 + 1;
       for ctr13_ovrtar = 1:size(param.numofSUMOontarget,2)
           
           if ctr13_ovrtar == 1
               tarIND = indE2;
           elseif ctr13_ovrtar == 2
               tarIND = indE1;
           elseif ctr13_ovrtar == 3
               tarIND =  param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1)) + (2*(param.numofSUMOontarget(2) + 1));
           else
               tarIND = tarIND + param.numofSUMOontarget(ctr13_ovrtar - 1) + 1;
           end
           
           %%% Reaction 13::  Pool->target  %target formation
           r(indrex13,1) = param.k13_e2e1tarform(1,ctr13_ovrtar); %%% r(13,1) = param.k(13,1);
           s(tarIND + 1,indrex13) = 1; %%% s(ntarget,13) = 1;
           indrex13 = indrex13 + 1 ; 
           
           %%% Reaction 14::  target -> pool %% Unsumoylated target degradation
           r(indrex13,1) = param.k14_e2e1tardegrade(1,ctr13_ovrtar)*y(tarIND + 1); %%% r(14,1) = param.k(14,1)*target;
           s(tarIND + 1,indrex13) = -1; %%% s(ntarget,14) = -1;
           indrex13 = indrex13 + 1 ;
           
           for ctr13_sumoONtar = 2:param.numofSUMOontarget(ctr13_ovrtar) + 1
               %%% Reaction 14::  target -> pool %% Sumoylated target degradation
               r(indrex13,1) = param.k14_e2e1tardegrade(ctr13_sumoONtar,ctr13_ovrtar)*y(tarIND + ctr13_sumoONtar); %%% r(14,1) = param.k(14,1)*target;
               s(tarIND + ctr13_sumoONtar ,indrex13) = -1; %%% s(ntarget,14) = -1;

               indrex13 = indrex13 + 1 ;  
           end
       end
         
         
 %%        
         
dy = s*r;

return;
