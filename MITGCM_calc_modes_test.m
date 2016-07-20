clear
base = '/import/c/w/jpender/MITgcm/simulations/Lighthill/out3_HLSfuncs/'
eval(['load ',base,'/matlab/data.mat']);
MODEL=MITGCM_calc_modes(MODEL,5,2*pi/(12.4*3600));
%% lets fit to a time series

ufile='/import/c/w/jpender/MITgcm/simulations/Lighthill/out3_HLSfuncs/netcdf/UVEL.nc';T=nc_varget(ufile,'T');
clear as*
for ii = 100
for jj = 40
 nWater = find(floor(MODEL.H(jj,ii) ./ MODEL.Z), 1, 'last');
 u=(nc_varget(ufile,'UVEL',[0,0,jj-1,ii-1],[-1,-1,1,1]));
 psi=sq(MODEL.psip(1:nWater,:,nWater));

 for tdx = 1:length(T) 
  a_p(tdx,:)=psi(:,:)'*u(tdx,1:nWater)'/nWater;
 end

end % ii
end % jj

figure(1);clf;mdx=2;
 subplot(2,1,1);plot(a_p(:,mdx)*psi(1,mdx),'k.-');hold on
                plot(u(:,1),'r');axis tight;
                legend(['upper level velocity from mode ',num2str(mdx)],...
                        'raw upper level velocity',3)
figure(2);clf;
mdx1=2;mdx2=3;
tmp=cumsum(psi(:,mdx1).*psi(:,mdx2).*MODEL.delZ(1:nWater));
    tmp=tmp/max(abs(tmp));
plot(tmp,'k.-');hold on;plot(tmp*0,'k--');text(nWater,.05,num2str(tmp(end)),'horizontalal','center')
title('modes are orthogonal if \int_0^H \psi_m \psi_n = 0')