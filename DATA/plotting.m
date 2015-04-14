myColors={'r','k','b','g','m','c','y'};
Pspan = [.0001 .001 .01 .1 1 10 100];
load('enp0001.dat')
load('enp001.dat')
load('enp01.dat')
load('enp1.dat')
load('en1.dat')
load('en10.dat')
load('en100.dat')

load('atm0001.dat')
load('atm001.dat')
load('atm01.dat')
load('atmp1.dat')
load('atm1.dat')
load('atm10.dat')
load('atm100.dat')

load('sp0001.dat')
load('sp001.dat')
load('sp01.dat')
load('sp1.dat')
load('s1.dat')
load('s10.dat')
load('s100.dat')

for o = 1:7
Pin = Pspan(o);

fOut=sprintf('T_%03d.mat',Pin);
load(fOut)
fOut=sprintf('Sstar_%03d.mat',Pin);
load(fOut)
fOut=sprintf('Estar_%03d.mat',Pin);
load(fOut)
fOut=sprintf('Hstar_%03d.mat',Pin);
load(fOut)
fOut=sprintf('Z_%03d.mat',Pin);
load(fOut)
fOut=sprintf('w_%03d.mat',Pin);
load(fOut)
fOut=sprintf('cstar_%03d.mat',Pin);
load(fOut)
fOut=sprintf('gamma_%03d.mat',Pin);
load(fOut)
fOut=sprintf('ecomp_%03d.mat',Pin);
load(fOut)
fOut=sprintf('astar_%03d.mat',Pin);
load(fOut)


% figure(1)
% hold on
% plot(T,ecomp,myColors{o})
% ylabel('Ratio of translational to electronic energy modes - e_{el}/e_t')
% xlabel('T - Temperature [K]')
% xlim([Tmin Tmax])
% ylim([0 .3])
% TEMPERATURE AND SPECIFIC TOTAL INTERNAL ENERGY
% figure(2)
% hold on
% plot(T,Z,'Color',myColors{o})
% ylabel('Z - Compressibility Factor')
% xlabel('T - Temperature [K]')
% xlim([Tmin Tmax])
% ylim([0 4])
% 
% TEMPERATURE AND DIMENSIONLESS SPECIFIC INTERNAL ENERGY
% figure(3)
% hold on
% plot(T,Estar,'Color',myColors{o})
% ylabel('E^* = (ZE_{tot})/(RT) - Dimensionless Specific Internal Energy')
% xlabel('T - Temperature [K]')
% xlim([Tmin Tmax])
% ylim([0 50])
% 
% 
% figure(4)
% hold on
% plot(T,astar,myColors{o})
% ylabel('a^* = a/a_{ideal} - Dimensionless speed of sound')
% xlabel('T - Temperature [K]')
% xlim([Tmin Tmax])
% 
% 
% % TEMPERATURE AND SPECIFIC TOTAL ENTHALPY
% figure(5)
% hold on
% plot(T,Sstar,'Color',myColors{o})
% ylabel('S^* = (ZS_{tot})/R - Dimensionless Specific Entropy')
% xlabel('T - Temperature [K]')
% xlim([Tmin Tmax])
% 
% figure(6)
% hold on
% plot(T,cstar,'Color',myColors{o})
% ylabel('c_v^* = (Zc_v)/R - Dimensionless specific heat ')
% xlabel('T - Temperature [K]')
% xlim([Tmin Tmax])
% 
figure(7)
hold on
plot(T,gamma,'Color',myColors{o})
ylabel('\gamma = c_p/c_v - Ratio of specific heats')
xlabel('T - Temperature [K]')
xlim([Tmin Tmax])

% TEMPERATURE AND SPECIFIC TOTAL ENTHALPY
% figure(8)
% hold on
% plot(T,hstar,myColors{o})
% ylabel('H^* = (ZH_{tot})/(RT) - Dimensionless specific enthalpy')
% xlabel('T - Temperature [K]')
% xlim([2000 15000])

% figure(9)
% hold on
% plot(T,RHOT,myColors{o})
% ylabel('\rho^* = \rho/\rho_0 - Dimensionless density')
% xlabel('T - Temperature [K]')
% xlim([2000 15000])
end
%% asdf
% figure(1)
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% figure(4)
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% figure(6)
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% figure(7)
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% figure(8)
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
figure(7)
legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% 
% figure(2)
% plot(atm0001(:,2),atm0001(:,3),'or')
% plot(atm001(:,2),atm001(:,3),'ok')
% plot(atm01(:,2),atm01(:,3),'ob')
% plot(atmp1(:,1),atmp1(:,2),'og')
% plot(atm1(:,1),atm1(:,2),'om')
% plot(atm10(:,1),atm10(:,2),'oc')
% plot(atm100(:,1),atm100(:,2),'oy')
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% figure(3)
% plot(enp0001(:,2),enp0001(:,3),'or')
% plot(enp001(:,2),enp001(:,3),'ok')
% plot(enp01(:,2),enp01(:,3),'ob')
% plot(en1(:,2),en1(:,3),'og')
% plot(enp1(:,2),enp1(:,3),'om')
% plot(en10(:,2),en10(:,3),'oc')
% plot(en100(:,2),en100(:,3),'oy')
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
% 
% figure(5)
% plot(sp0001(:,2),sp0001(:,3),'or')
% plot(sp001(:,2),sp001(:,3),'ok')
% plot(sp01(:,2),sp01(:,3),'ob')
% plot(sp1(:,2),sp1(:,3),'og')
% plot(s1(:,2),s1(:,3),'om')
% plot(s10(:,2),s10(:,3),'oc')
% plot(s100(:,2),s100(:,3),'oy')
% legend('P = 0.0001 atm','P = 0.001 atm','P = 0.01 atm','P = 0.1 atm','P = 1 atm','P = 10 atm','P = 100 atm')
