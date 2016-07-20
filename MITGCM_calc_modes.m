function MODEL=MITGCM_calc_modes(MODEL,nmodes,omega)
%%

rho=(sw_pden(34+0*real(MODEL.OB.ET(1,:,1)),real(MODEL.OB.ET(1,:,1)),MODEL.Z',30));

Nx=MODEL.Nx;  Ny=MODEL.Ny;  Nz=MODEL.Nz;Hmax=MODEL.H_max;z=MODEL.Z;

%% Construct the basis functions. 

N2 = MODEL.N2;%
z  = MODEL.Z; 
 dz=MODEL.delZ;
 P = sw_pres(z,MODEL.reflat);

psip = nan*ones([MODEL.Nz  nmodes  MODEL.Nz ]);

for nzi = 4:Nz ;
[psidw,psidp,ce]       = dynmodes_hls(N2(1:nzi),P(1:nzi),sw_f(MODEL.reflat),2*pi/(12.4*3600));
psidp=[(1+0*psidp(2:end,1)) psidp(2:end,:)];psidw=[(1+0*psidw(2:end,1)) psidw(2:end,:)];% add BT mode
ce=[sqrt(9.8*z(nzi));ce];
% normalize pmodes
    A=repmat(nansum(psidp.^2.*dz(1),1)./z(nzi),[nzi 1]).^(1/2);
    A(A==0)=Inf;
    psidp=psidp./A;
    psidp(:,psidp(end,:)<0)=-psidp(:,psidp(end,:)<0);

%populate psi matrix
psip(1:nzi,1:min(nzi,nmodes),nzi)=psidp(:,1:min(nzi,nmodes));
end
MODEL.psip=psip;
%
