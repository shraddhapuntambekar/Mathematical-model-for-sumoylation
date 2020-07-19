function y0= initial(param,y0Targets)

    getname;

    y0 = zeros(1, param.n_species);
    
    y0(1,npresumo)=0;
    y0(1,nsumo)= 10;
    y0(1,nsenp)= 1;
    y0(1,ne3)= 1;
    
    y0(1,param.numOFconstspecies + 1) = y0Targets(1);  %% E20
    y0(1,param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1)) + 1) = y0Targets(2);  %% E10

    tempindtar =  param.numOFconstspecies + (2*(param.numofSUMOontarget(1) + 1)) + (2*(param.numofSUMOontarget(2) + 1));
    for indTAR = 3:size(param.numofSUMOontarget,2) 
        y0(tempindtar + 1) = y0Targets(indTAR);
        
        tempindtar = tempindtar +  (param.numofSUMOontarget(indTAR) + 1);
    end
    
    y0=y0/param.cref; 

return;
