function [wmodes,pmodes,ce]=dynmodes(Nsq,p,f,omega)
% DYNMODES calculates ocean dynamic vertical modes
%  taking a column vector of Brunt-Vaisala values (Nsq) at
%  different pressures (p) and calculating some number of 
%  dynamic modes (nmodes). 
%  Note: The input pressures need not be uniformly spaced, 
%    and the deepest pressure is assumed to be the bottom.
%
%  USAGE: [wmodes,pmodes,ce]=dynmodes(Nsq,p,nmodes);
%                               or
%                            dynmodes;  % to show demo 
%
%     Inputs: 	Nsq = column vector of Brunt-Vaisala buoyancy frequency (s^-2)
%		    	  p = column vector of pressure (decibars)
%           nmodes = number of vertical modes to calculate 
%  
%       Outputs:   wmodes = vertical velocity structure
%                  pmodes = horizontal velocity structure
%                      ce = modal speed (m/s)
%  developed by J. Klinck. July, 1999
%  send comments and corrections to klinck@ccpo.odu.edu
%  modifiled by H. Simmons 2012
rho0=1028;

%    convert to column vector if necessary
[m,n] = size(p);
if n == 1
   p=p';
end
[m,n] = size(Nsq);
if n == 1
   Nsq=Nsq';
   n=m;
end

%                 check for surface value
if p(1) > 0
%             add surface pressure with top Nsq value
    z(1)=0;
    z(2:n+1)=-p(1:n);
    N2(1)=Nsq(1);
    N2(2:n+1)=Nsq(1:n);
    nz=n+1;
else
    z=-p;
    N2=Nsq;
    nz=n;
end

%          calculate depths and spacing
%        spacing
dz(1:nz-1)=z(1:nz-1)-z(2:nz);
%        midpoint depth
zm(1:nz-1)=z(1:nz-1)-.5*dz(1:nz-1)'';
%        midpoint spacing
dzm=zeros(1,nz);
dzm(2:nz-1)=zm(1:nz-2)-zm(2:nz-1);
dzm(1)=dzm(2);
dzm(nz)=dzm(nz-1);

%        get dynamic modes
D2 = zeros(nz,nz);
B = zeros(nz,nz);
%             create second derivative matrix   
for i=2:nz-1
  D2(i,i-1) = -1/(dz(i-1)*dzm(i));
  D2(i,i  ) = +1/(dz(i-1)*dzm(i))  + 1/(dz(i)*dzm(i));
  D2(i,i+1) = -1/(dz(i)*dzm(i));
end

for i=1:nz
  % B(i,i)=N2(i); % orig
  B(i,i)=N2(i);%/(omega^2-f^2);
end

%             set boundary conditions
D2(1 ,1) =-1.; % rigid lid   % HLS
D2(nz,1) =-1.; % flat bottom % HLS
%
[wmodes,e] = eig(D2,B);

%          extract eigenvalues
e=diag(e);
%
ind=find(imag(e)==0);
e=e(ind);
wmodes=wmodes(:,ind);
%
ind=find(e>=1.e-10);
e=e(ind);
wmodes=wmodes(:,ind);
%
[e,ind]=sort(e);
wmodes=wmodes(:,ind);

nm=length(e);
ce=1./sqrt(e);
%                   create pressure structure
pmodes=zeros(size(wmodes));

for i=1:nm
%           calculate first deriv of vertical modes
  pr=diff(wmodes(:,i));   
  pr(1:nz-1)= pr(1:nz-1)./dz(1:nz-1)';
  pr=pr*rho0*ce(i)*ce(i);
%       linearly interpolate back to original depths
  pmodes(2:nz-1,i)=.5*(pr(2:nz-1)+pr(1:nz-2));
  pmodes(1,i)=pr(1);
  pmodes(nz,i)=pr(nz-1);
end


