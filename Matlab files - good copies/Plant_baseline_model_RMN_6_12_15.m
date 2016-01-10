%MINIMAL PLANT MODEL - BASED ON NISBET NOTES OF 6/12/15

%CODE LAST UPDATED 5/25/15

clear all

%CHOICE OF IMPLEMENTATIONS OF SYNTHESIZING UNIT - COMMENT THE ONE NOT USED
synth = @(x,y) min(x,y);      % MINIMUM RULE
%synth = @(x,y) x*y*(x+y)/(x^2+ x*y + y^2);    %PARALLEL, COMPLEMENTARY SU


%Run specification
dt=0.0001;		% time step.
tmax=10;     % run length
t=0:dt:tmax;  %time vector.

%model parameters
alphaC=1;  alphaN=1; gammaR=0.1; gammaS=0.1; betaS=0.25; betaR=1.0; sigmaR=0.0; sigmaS=0.0;


%initialization - somewhat arbitrary
S(1)=0.2;  R(1)=1.0; 
UC(1) = alphaC*S(1);  UN(1) = alphaN*R(1); 
rhoC(1)=UC(1); rhoN(1)=UN(1); 
TS(1) = gammaS*S(1);  TR(1) = gammaR*R(1);
rS(1) = sigmaS*betaS*gammaS*S(1); rR(1)=sigmaR*betaR*gammaR*R(1);
QS(1) = 0;  QR(1)=0;
wN(1) = 0;
wC(1) = 0;  

%updating
for i=2:length(t),	
   UC(i)=alphaC*S(i-1);    UN(i)=alphaN*R(i-1);
   TS(i)=gammaS*S(i-1);    TR(i)=gammaR*R(i-1);
   rS(i)=sigmaS*betaS*gammaS*S(i-1); rR(i)=sigmaR*betaR*gammaR*R(i-1);
   QS(i)=synth(UC(i), (rhoN(i-1)+rS(i))/betaS);
   QR(i)=synth(rhoC(i-1),(UN(i)+rR(i))/betaR);
   rhoC(i)=max(UC(i)-QS(i),0); 
   rhoN(i)=max(UN(i)+rR(i)-betaR*QR(i),0);
   wC(i)=max(rhoC(i)-QR(i),0);
   wN(i)=max(rhoN(i)+rS(i)*betaS-QS(i)*betaS,0);
   S(i)=S(i-1)+dt*(QS(i)-TS(i));
   R(i)=R(i-1)+dt*(QR(i)-TR(i));
   end


%Estimate final slope
M=length(t);
R(M);
S(M);
slope=(log(R(M))-log(R(M-10)))/(10*dt)

%Plots.
subplot(2,2,1)
plot(t,log(R),t, log(S),'r')
%hold on
%plot(t,log(S), 'r')
title('Root (blue) and shoot (red) biomass - log plot');
xlabel('Time');
ylabel('log(S) and log(H)');
%hold off

subplot(2,2,2)
plot(t,rhoN./UN)
hold on
plot(t,rhoC./UC,'r')
hold off
title('Rejection flux to uptake ratio.  red: carbon, blue: nitrogen');
xlabel('time')
ylabel('rhoN and rhoC')

subplot(2,2,3)
plot(t,wN./UN)
hold on
plot(t,wC./UC,'r')
hold off
title('waste fluxe to uptake ratio - red: carbon, blue: nitrogen');
xlabel('time')
ylabel('wN and wC')

