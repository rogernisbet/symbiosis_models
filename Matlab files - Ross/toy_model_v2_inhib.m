% MINIMAL CORAL MODEL
% BASED ON MINIMAL PLANT MODEL - NISBET NOTES OF 5/25/15

% CODE LAST UPDATED 9/18/15 - RC

clear all

%CHOICE OF IMPLEMENTATIONS OF SYNTHESIZING UNIT - COMMENT THE ONE NOT USED
%synth = @(x,y) min(x,y);      % MINIMUM RULE
synth = @(x,y) x*y*(x+y)/(x^2 + x*y + y^2);    %PARALLEL, COMPLEMENTARY SU


%Run specification
dt=0.001;		% time step.
tmax=50;     % run length
t=0:dt:tmax;  %time vector.

% Define environmental vectors (Light (L), Prey (X), DIN (N))
L=linspace(8e7,8e7,length(t)); % constant light vector
%L=linspace(8e7,8e12,length(t)); % increasing light vector
%L=100*sin(pi*t) + 100; % sinusoidal light vector
%L=[linspace(8e13,8e8,length(t)/2),linspace(8e7,8e7, length(t)/2+1)]; % decreasing to zero light vector
%L=[linspace(8e9,8e9,length(t)/4),linspace(8e13,8e13,length(t)/4),linspace(5e6,5e6,length(t)/4),linspace(0,0,length(t)/4+1)];
X=linspace(1,1,length(t)); % food vector  % less food --> higher S/H
N=linspace(0.01,0.01,length(t)); % DIN vector

%model parameters
gammaR = 0.2; % biomass loss rate of coral per C-mol
gammaS = 0.2; % biomass loss rate of symbiont per C-mol
nnS = 0.15; % N:C ratio in symbiont
nnR = 0.2; % N:C ratio in coral
nnX = 0.2; % N:C ratio in prey
yCL = 0.05; % yield of carbon from light (photons)
Cmax= 6; % maximum carbon production per C-mol of symbiont
ds=1/(5e-6); % symbiont C-mol / cross sectional area
jIp=50e6;  % inhibition scaling parameter -- what is this? photons / sec?
Ipswitch=1; % photoinhibition on/off switch

%initialization values
S(1)=0.01;  R(1)=1; % specify symbiont and coral starting biomass
UL(1) = S(1) / ds * L(1);  % uptake of light prop. to sym biomass and external light
UX(1) = R(1) * X(1);
UN(1) = R(1) * N(1);
UC(1) = 0 ;
QS(1) = 0;  
QR(1) = 0;
rhoC(1) = UC(1); 
rhoN(1) = UN(1); 
TS(1) = gammaS*S(1);  
TR(1) = gammaR*R(1);
 

%updating
for i = 2:length(t),	
   %UL(i) = S(i-1) * L(i);  % uptake of light (symbiont)
   %UC(i) = synth(UL(i) * yCL, TR(i-1));  % uptake of carbon (symbiont)
   UL(i) = L(i) / ds * S(i-1);  % total uptake of light at time t
   JIp(i) = jIp * S(i-1); % inhibition scaling -- should be multiplied by S biomass?
   B1(i) = 1/TR(i-1);  % 1/total flux of carbon to symbiont photosynthesis
   C1(i) = 1/(yCL * UL(i));  % 1/total flux of light to symbiont photosynthesis
   D1(i) = 1/(TR(i-1)+yCL*UL(i));
   inhibp(i) = (1+Ipswitch*UL(i)/JIp(i));
   UC(i) = 1 / ((1/(Cmax*S(i-1))*inhibp(i))+(B1(i)+C1(i)-D1(i)));
   UX(i) = R(i-1) * X(i);  % uptake of food (coral)
   UN(i) = R(i-1) * N(i);  % uptake of DIN (coral)
   TS(i) = gammaS * S(i-1);  % turnover of symbiont
   TR(i) = gammaR * R(i-1);  % turnover of coral
   QS(i) = synth(UC(i), (rhoN(i-1))/nnS);  % formation of symbiont biomass
   QR(i) = synth(rhoC(i-1) + UX(i-1), (UN(i) + nnX * UX(i-1))/nnR);  % formation of coral biomass
   rhoC(i) = max(UC(i) - QS(i), 0);  % rejection of carbon from symbiont to coral
   rhoN(i) = max(UN(i) + nnX * UX(i) - nnR * QR(i), 0);  % rejection of nitrogen from coral to symbiont
   S(i) = S(i-1) + dt * (QS(i) - TS(i));  % change in symbiont biomass
   R(i) = R(i-1) + dt * (QR(i) - TR(i));  % change in coral biomass
end


%Estimate final slope (growth rate)
M=length(t);
R(M);
S(M);
slope=(log(R(M))-log(R(M-10)))/(10*dt)

%Plots.
subplot(5,1,1)
plot(t, log(R), t, log(S), 'r')
%hold on
%plot(t,log(S), 'r')
title('Coral (blue) and symbiont (red) biomass - log plot');
xlabel('Time');
ylabel('log(S) and log(R)');
%hold off

subplot(5,1,2)
plot(t, S./R, 'black')
title('Symbiont to host ratio');
xlabel('Time')
ylabel('S/R')

subplot(5,1,3)
plot(t, log(L+1), 'black')
title('Light');
xlabel('Time')
ylabel('Light')

subplot(5,1,4)
plot(t, N, 'black')
title('DIN');

subplot(5,1,5)
plot(t, X, 'black')
title('food');

%subplot(3,2,3)
%plot(t,rhoN./UN)
%hold on
%plot(t,rhoC./UC,'r')
%hold off
%title('Rejection flux to uptake ratio.  red: carbon, blue: nitrogen');
%xlabel('time')
%ylabel('rhoN and rhoC')

%subplot(3,2,4)
%plot(t,wN./UN)
%hold on
%plot(t,wC./UC,'r')
%hold off
%title('waste fluxe to uptake ratio - red: carbon, blue: nitrogen');
%xlabel('time')
%ylabel('wN and wC')

