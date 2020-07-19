
n_species =18; 
dum1=num2cell(zeros(n_species ,1));
for n = 1:n_species
    dum1{n} = n; %value assignment to each cell.. dum{1}=1; dum{2}=2...
end

[npresumo,nsumo,nsenp,ne10,nse10,ne20,ne21,nse20,nse21,nt10,nt11,nt20,nt21,npresumosenp,nt11senp,nt21senp,ne21senp,nsumoe21senp]= deal(dum1{:});

clear  n dum1