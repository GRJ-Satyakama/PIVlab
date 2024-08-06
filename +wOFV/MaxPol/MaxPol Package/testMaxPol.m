clear variables
close all
clc

%Test MaxPol compared to some naive gradient methods, and play with its
%parameters for getting image gradients. The test image is a 2D sinusoid so
%that the gradients have analytical solutions.

N=512;
[testIm,trueGd]=sinusoidal_gratting_2D(.4,0,N);
lfix=3;
Pfix=2;

q=1;
figure(q)
imshow(testIm,[min(testIm(:)),max(testIm(:))])
q=q+1;

trueGd_c=cat(3,trueGd.centralized{2},trueGd.centralized{1})*2/(N-1);
trueGd_s=cat(3,trueGd.staggered{2},trueGd.staggered{1})*2/(N-1);
minG=min(min(trueGd_c(:,:,1)));
maxG=max(max(trueGd_c(:,:,1)));

% figure(q)
% imshow(trueGd_s(:,:,1),[])
% q=q+1;
% figure(q)
% imshow(trueGd_s(:,:,2),[])
% q=q+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start with noise-free image

sparseflag=1;
symflag='true';
nod='centralized';

disp('Naive OFgrad function (central difference):')
tic
testGd=OFgrad(testIm);
toc
OFerr=RMSEvec(testGd,trueGd_c);
disp(['Error: ' num2str(OFerr)])
disp(' ')
% figure(q)
% subplot(1,2,1)
% imshow(trueGd_c(:,:,1),[minG,maxG])
% subplot(1,2,2)
% imshow(testGd(:,:,1),[minG,maxG])
% q=q+1;

llist=1:8;
errvec=zeros(length(llist),1);
c=1;
for l=llist
    disp(['MaxPol, l=' num2str(l) ':'])
    D=derivmtx(l,Pfix,1,N,nod,sparseflag,symflag);
    tic
    testGd=cat(3,D*testIm,testIm*D');
    toc
    disp(['Error: ' num2str(RMSEvec(testGd,trueGd_c))])
    disp(' ')
%     figure(q)
%     subplot(1,2,1)
%     imshow(trueGd_c(:,:,1),[minG,maxG])
%     subplot(1,2,2)
%     imshow(testGd(:,:,1),[minG,maxG])
%     q=q+1;
    errvec(c)=RMSEvec(testGd,trueGd_c);
    
    c=c+1;
end
figure(q)
semilogy(llist,errvec)
hold on
semilogy([llist(1),llist(end)],[OFerr,OFerr],'--')
q=q+1;
% 
% Plist=(1:lfix)*2;
% errvec=zeros(length(Plist),1);
% c=1;
% for P=Plist
%     disp(['MaxPol, P=' num2str(P) ':'])
%     D=derivmtx(lfix,P,1,N,nod,sparseflag,symflag);
%     tic
%     testGd=cat(3,D*testIm,testIm*D');
%     toc
%     disp(['Error: ' num2str(RMSEvec(testGd,trueGd_c))])
%     disp(' ')
% %     figure(q)
% %     subplot(1,2,1)
% %     imshow(trueGd_c(:,:,1),[minG,maxG])
% %     subplot(1,2,2)
% %     imshow(testGd(:,:,1),[minG,maxG])
% %     q=q+1;
%     errvec(c)=RMSEvec(testGd,trueGd_c);
%     
%     c=c+1;
% end
% figure(q)
% semilogy(Plist,errvec)
% hold on
% semilogy([Plist(1),Plist(end)],[OFerr,OFerr],'--')
% q=q+1;
% 
% figure(q)
% imshow(log(abs(trueGd_c(:,:,1)-testGd(:,:,1))),[])
% colormap('parula')
% q=q+1;

%Test with Tom's differentiation settings
D=derivmtx(3,2,1,N,nod,sparseflag,symflag);
S=derivmtx(2,2,0,N,nod,sparseflag,symflag);

disp('MaxPol, Tom''s settings:')
tic
testGd=cat(3,D*testIm*S,S'*testIm*D');
toc
disp(['Error: ' num2str(RMSEvec(testGd,trueGd_c))])
disp(' ')
semilogy([llist(1),llist(end)],[RMSEvec(testGd,trueGd_c),RMSEvec(...
    testGd,trueGd_c)],'--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add noise
disp('*************************************************************')
disp('Noisy image')
disp('*************************************************************')

nlev=.02*max(testIm(:));
testIm=testIm+nlev*randn(N);
figure(q)
imshow(testIm,[min(testIm(:)),max(testIm(:))])
q=q+1;

disp('Naive OFgrad function (central difference):')
tic
testGd=OFgrad(testIm);
toc
OFerr=RMSEvec(testGd,trueGd_c);
disp(['Error: ' num2str(OFerr)])
disp(' ')
% figure(q)
% subplot(1,2,1)
% imshow(trueGd_c(:,:,1),[minG,maxG])
% subplot(1,2,2)
% imshow(testGd(:,:,1),[minG,maxG])
% q=q+1;

llist=1:8;
errvec=zeros(length(llist),1);
c=1;
for l=llist
    disp(['MaxPol, l=' num2str(l) ':'])
    D=derivmtx(l,Pfix,1,N,nod,sparseflag,symflag);
    tic
    testGd=cat(3,D*testIm,testIm*D');
    toc
    disp(['Error: ' num2str(RMSEvec(testGd,trueGd_c))])
    disp(' ')
%     figure(q)
%     subplot(1,2,1)
%     imshow(trueGd_c(:,:,1),[minG,maxG])
%     subplot(1,2,2)
%     imshow(testGd(:,:,1),[minG,maxG])
%     q=q+1;
    errvec(c)=RMSEvec(testGd,trueGd_c);
    
    c=c+1;
end
figure(q)
semilogy(llist,errvec)
hold on
semilogy([llist(1),llist(end)],[OFerr,OFerr],'--')
q=q+1;
% 
% errvec=zeros(length(Plist),1);
% c=1;
% for P=Plist
%     disp(['MaxPol, P=' num2str(P) ':'])
%     D=derivmtx(lfix,P,1,N,nod,sparseflag,symflag);
%     tic
%     testGd=cat(3,D*testIm,testIm*D');
%     toc
%     disp(['Error: ' num2str(RMSEvec(testGd,trueGd_c))])
%     disp(' ')
% %     figure(q)
% %     subplot(1,2,1)
% %     imshow(trueGd_c(:,:,1),[minG,maxG])
% %     subplot(1,2,2)
% %     imshow(testGd(:,:,1),[minG,maxG])
% %     q=q+1;
%     errvec(c)=RMSEvec(testGd,trueGd_c);
%     
%     c=c+1;
% end
% figure(q)
% semilogy(Plist,errvec)
% hold on
% semilogy([Plist(1),Plist(end)],[OFerr,OFerr],'--')
% q=q+1;

% figure(q)
% imshow(log(abs(trueGd_c(:,:,1)-testGd(:,:,1))),[])
% colormap('parula')
% q=q+1;

%Test with Tom's differentiation settings
tic
D=derivmtx(3,2,1,N,nod,sparseflag,symflag);
S=derivmtx(2,2,0,N,nod,sparseflag,symflag);
toc

disp('MaxPol, Tom''s settings:')
tic
testGd=cat(3,D*testIm*S,S'*testIm*D');
toc
disp(['Error: ' num2str(RMSEvec(testGd,trueGd_c))])
disp(' ')
semilogy([llist(1),llist(end)],[RMSEvec(testGd,trueGd_c),RMSEvec(...
    testGd,trueGd_c)],'--')