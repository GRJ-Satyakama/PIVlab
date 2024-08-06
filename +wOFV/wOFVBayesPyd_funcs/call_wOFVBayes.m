function [valF,valG]=call_wOFVBayes(Theta,I0,I1,scale,eta,vartheta,Fw,FwInv,Ni,...
    regscheme,Dmat,mask)

%Creates persistent variables for wavelet optical flow function call 
%(wOFVBayes_min)

persistent I0S I1S scaleS FwS FwInvS NiS regschemeS DmatS maskS etaS varthetaS

if nargin>=2
    I0S=I0;
    I1S=I1;
    scaleS=scale;
    FwS=Fw;
    FwInvS=FwInv;
    NiS=Ni;
    regschemeS=regscheme;
    DmatS=Dmat;
    maskS=mask;
    etaS = eta;
    varthetaS = vartheta;
    return
end

[valF,valG]=wOFVBayes_min(Theta,I0S,I1S,scaleS,etaS,varthetaS,FwS,FwInvS,...
    NiS,regschemeS,DmatS,maskS);