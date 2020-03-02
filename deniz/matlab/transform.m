clear all;
close all;
load('FilterCoefficients.mat');
fileID = fopen('myOutFile.txt','r');
ccc=fscanf(fileID,'%f');
fileID = fopen('freqs2.txt','w');
fileID2 = fopen('J05.txt','w');
eta=0;
c=0;
tau=0;
sigmaInf=1;
mu_0=1.2566370614359173e-06;
r  = sqrt( 200^2. + 0^2. + 0^2.);
x=200;
time=logspace(-3,1,10);
for tt=1:length(time)
    t=time(tt);
    theta=sqrt(mu_0*sigmaInf/(4*t));
    errfunc=erfc(theta*r);
    firstterm=( ( (4/sqrt(pi))*(theta^3)*(r^3) + ((6/sqrt(pi))*theta*r) )*exp(-theta^2*r^2)+3*errfunc );
    last=((((4/sqrt(pi))*(theta^3)*(r^3)+(2/sqrt(pi))*theta*r))*exp(-theta^2*r^2)+errfunc);
    ee(tt)=(1/(4*pi*sigmaInf*r^3))*(firstterm*(x^2/r^2)-last);
    for neww=1:N_J05
        w=(root_J05)^(neww-N0_J05)/t;
        w=w*(9.054710985091e-01);
        freqs(neww,tt)=w;
        sigma=sigmaInf*(1 - eta/(1 + (1i*2*pi*w*tau)^c));
        k  = sqrt(-1i*w*mu_0*sigma);
        front = 1 * 1 / (4. * pi * sigma * r^3) * exp(-1i*k*r);
        mid   = -k^2 * r^2 + 3*1i*k*r + 3;
        WH240(neww)=(1/(1i*w))*front*((x^2 / r^2)*mid + (k^2 * r^2 -1i*k*r-1));
        %WH240(neww)=front*((x^2 / r^2)*mid + (k^2 * r^2 -1i*k*r-1));
    end
    Efield(tt)=sqrt(2/pi)*imag(WH240)*J05/t;
    
end
fprintf(fileID,'%.8g\n',freqs(:,[10]));
fclose(fileID);
fprintf(fileID2,'%.8g\n',J05);
fclose(fileID2);
semilogx(time,abs(Efield),'--','LineWidth',5,'DisplayName','W&H 2.40')
hold on
semilogx(time,ee,'DisplayName','W&H 2.50')
xlabel('time')
ylabel('E')
legend
hold off
