function y0= getinitialCLOSED(param)

getnames;

y0(1,npresumo) =1;
y0(1,nsumo)= 0;
y0(1,nsenp)= 1;
y0(1,ne10)=1;
y0(1,nse10)=0;
y0(1,ne20)=1;
y0(1,ne21)=0;
y0(1,nse20)=0;
y0(1,nse21)=0;
y0(1,nt10)=1;
y0(1,nt11)=0;
y0(1,nt20)=0;
y0(1,nt21)=0;
y0(1,npresumosenp)=0;
y0(1,nt11senp)=0;
y0(1,nt21senp)=0;
y0(1,ne21senp)=0;
y0(1,nsumoe21senp)=0;

y0=y0/param.cref; 

return;