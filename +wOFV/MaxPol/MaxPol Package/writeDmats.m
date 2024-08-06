clear variables
close all
clc

%Write D-matrices for given image sizes
destdir='../../Optical Flow Velocimetry/WOF_mk3/Filter matrices/Diff/';

lfix=1;
Pfix=2;

sparseflag=1;
symflag='true';
nod='centralized';

dord=1;

for N=2.^(5:11)
    mkdir([destdir num2str(N)])
    
    Dmat=derivmtx(lfix,Pfix,dord,N,nod,sparseflag,symflag);
    
    save([destdir num2str(N) '/Dmat.mat'],'Dmat')
end