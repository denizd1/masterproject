
function [V,Ephi,Bz,Br] = Layer(sigma1,sigma2,sigma3,hTX,hRX,thickness,a,r,I,t,J);

load('FilterCoefficients.mat');

[Dgs] = GaverStehfest(J);

[sigmaBelow,dBelow,sigmaAbove,dAbove,NRX]=dummylayering(sigma1,sigma2,sigma3,hTX,hRX,thickness);

N=length(sigmaBelow);
M=length(sigmaAbove);

epr=10;

    for h=1:N_J1;
        l=(root_J1)^(h-N0_J1)/a;
               
%         for neww=1:length(Dgs);
%             s=neww*((log(2))/t);
 
        for neww=1:N_J05;
            w=(root_J05)^(neww-N0_J05)/t;
            w=w*(9.054710985091e-01);
%            s=-i*w;

%             thetaQ(N)=sqrt(l^2+s*mu0*sigmaBelow(N)+s^2*mu0*ep0*epr);
%             Q(N)=mu0/thetaQ(N);
% 
%             thetaP(M)=sqrt(l^2+s*mu0*sigmaAbove(M)+s^2*mu0*ep0*epr);
%             P(M)=mu0/thetaP(M);
            
            
 %Water Layer Impedances
            
%             for trev=1:N-1;
%                 tt=N-trev;
%                 thetaQ(tt)=sqrt(l^2+s*mu0*sigmaBelow(tt)+s^2*mu0*ep0*epr);
%                 argQ(tt)=thetaQ(tt)*dBelow(tt);
%                 Q(tt)=(mu0/thetaQ(tt))*((thetaQ(tt)*Q(tt+1)+mu0*tanh(argQ(tt)))/(mu0+thetaQ(tt)*Q(tt+1)*tanh(argQ(tt))));
%             end;
%             
% 
%             for trev=1:M-1;
%                 tt=M-trev;
%                 thetaP(tt)=sqrt(l^2+s*mu0*sigmaAbove(tt)+s^2*mu0*ep0*epr);
%                 argP(tt)=thetaP(tt)*dAbove(tt);
%                 P(tt)=(mu0/thetaP(tt))*((thetaP(tt)*P(tt+1)+mu0*tanh(argP(tt)))/(mu0+thetaP(tt)*P(tt+1)*tanh(argP(tt))));
%             end;
%             
%             if hRX>hTX;
%                for lll=1:NRX-1;
%                 
%                 L(lll)=[Q(lll+1)*thetaQ(lll)*sech(dBelow(lll)*thetaQ(lll))]/[Q(lll+1)*thetaQ(lll)+mu0*tanh(dBelow(lll)*thetaQ(lll))];
%                 
%                end;
%             end;
%             
%             if hRX<hTX;
%                 for lll=1:NRX-1;
%                     
%                 L(lll)=[P(lll+1)*thetaP(lll)*sech(dAbove(lll)*thetaP(lll))]/[P(lll+1)*thetaP(lll)+mu0*tanh(dAbove(lll)*thetaP(lll))];
%                 
%                 end;
%             end;
%             
            
            theta=sqrt(l^2+s*mu0*sigma1+s^2*mu0*ep0*epr);
            
            V_ls(neww)=[2*pi]*[a^2]*[I/s]*[l]*[mu0/(2*theta)]*besselj(1,l*a); 
%             Ephi_ls(neww)=[a]*[I]*[l]*[(Q(1)*P(1))/(Q(1)+P(1))]*besselj(1,l*r)*prod(L);  
%             Bz_ls(neww)=[a]*[I]*[1/s]*[l^2]*[(Q(1)*P(1))/(Q(1)+P(1))]*besselj(0,l*r)*prod(L);
%             
%             if hRX>hTX;
%             Br_ls(neww)=-[a]*[I]*[mu0]*[1/s]*[l]*[(Q(1)*P(1))/(Q(1)+P(1))]*besselj(1,l*r)*prod(L)/Q(NRX); 
%             end;
%             
%             if hRX<hTX;
%             Br_ls(neww)=+[a]*[I]*[mu0]*[1/s]*[l]*[(Q(1)*P(1))/(Q(1)+P(1))]*besselj(1,l*r)*prod(L)/P(NRX); 
%             end;
%             
            
        end;
%             
%            V_lt(h)=V_ls*Dgs*((log(2))/t);
%             Ephi_lt(h)=Ephi_ls*Dgs*((log(2))/t);
%             Bz_lt(h)=Bz_ls*Dgs*((log(2))/t);
%             Br_lt(h)=Br_ls*Dgs*((log(2))/t);

            V_lt(h)=[sqrt(2/pi)*imag(V_ls)]*J05/t; 
%             Ephi_lt(h)=[sqrt(2/pi)*imag(Ephi_ls)]*J05/t;
%             Bz_lt(h)=[sqrt(2/pi)*imag(Bz_ls)]*J05/t;
%             Br_lt(h)=[sqrt(2/pi)*imag(Br_ls)]*J05/t;
            



    end;
    
    V=V_lt*J1/a;
%     Ephi=Ephi_lt*J1/a;
%     Bz=Bz_lt*J1/a;
%     Br=Br_lt*J1/a;
    
Ephi=0;
Bz=0;
Br=0;
    