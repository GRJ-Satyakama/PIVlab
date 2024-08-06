function err=RMSEvec(v1,v2)
%Compute root mean square error from two vector fields v1 and v2
difv=v1-v2;
normsq=difv(:,:,1).^2+difv(:,:,2).^2;

% size(normsq)
% length(normsq(:))

err=sqrt(sum(normsq(:))/length(normsq(:)));