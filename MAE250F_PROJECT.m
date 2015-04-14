% MAE 250F PROJECT
% CONSTANT PRESSURE
% Jon Van Lew, 2012

clear all
% close all
format short;format compact
clc

% INPUTS
Pspan = 1;%[.0001 .001 .01 .1 1 10 100];
%atm
Tmax = 15000; 							%K
Tmin = 2000; 							%K
N = 1000;



%% INITIALIZE ALL VARIABLES
% Species numbers in program
% 1 = N2
% 2 = N
% 3 = O2
% 4 = O
% 5 = NO
% 6 = e-
% 7 = O+
% 8 = N+
% 9 = NO+
Y = 8; % NUMBER OF SPECIES
% PROGRAMMING VALUES
z = (Tmax - Tmin)/(N-1);

% UNIVERSAL CONSTANTS
RGAS = 8314; 							%J/kmol-K
H = 6.625E-34; 							%M2-KG/S
K = 1.381E-23; 							%M2-KG/S2
NA = 6.022E23*1000; 					%/KMOL
% GAS PROPERTIES
mM = [28.0134,14.0067,31.9988,15.9994,30.0061,9.10938E-31*NA,15.9994,14.0067]; %KG/KMOL -- N2, N, O2, O, NO
M = mM/NA; 								%KG/MOLECULE
R = RGAS./mM;                           %J/KG-K
R0 = RGAS/(28.97);      %J/KG-K
% CHARACTERISTIC TEMPERATURES
% ROTATIONAL ENERGY
TH_R = [2.9,0,2.1,0,2.5]; 			% -- N2, N, O2, O, NO
% VIBRATIONAL ENERGY
TH_V = [3390,0,2270,0,2740]; 		% -- N2, N, O2, O, NO
% HEAT OF DISSOCATION
TH_D = [113000,0,59500,0,75500]; 	% -- N2, N, O2, O, NO
% IONIZATION ENERGY
TH_I = [181,168.8,142,158,108]*10^3; 	% -- N2, N, O2, O, NO
% STOICHIOMETRIC COEFFICIENTS
SIGMA = [2,0,2,0,1]; 					% -- N2, N, O2, O, NO


g(1,:) = [1,4,3,5,2,2,4,1];
g(2,:) = [0,10,2,3,2,0,10,3];
g(3,:) = [0,6,2,1,0,0,6,5];
g(4,:) = [0,0,0,5,0,0,0,5];
g(5,:) = [0,0,0,1,0,0,0,1];
g(6,:) = [0,0,0,0,0,0,0,5];
th_el(1,:) = [0,0,0,0,0,0,0,0];
th_el(2,:) = [0,27700,11390,228,174,0,38600,70.6];
th_el(3,:) = [0,41500,18990,326,0,0,58200,188.9];
th_el(4,:) = [0,0,0,22800,0,0,0,22000];
th_el(5,:) = [0,0,0,48600,0,0,0,47000];
th_el(6,:) = [0,0,0,0,0,0,0,67900];

cv_t=   zeros(1,Y);
cv_r=   zeros(1,Y);
cv_v=   zeros(1,Y);
cv_el=  zeros(1,Y);
cv=     zeros(N,Y);
gamma=  zeros(1,Y);
RT=     zeros(N,1);
Q=      zeros(N,Y);
T=      zeros(N,1);
Q_t=    zeros(1,Y);
Q_r=    zeros(1,Y);
Q_v=    zeros(1,Y);
Q_el=   zeros(1,Y);
KP=     zeros(N,Y);
PARTP=  zeros(N,Y);
P_TOT=  zeros(N,1);
X=      zeros(N,Y);
RHO=    zeros(N,Y);
RHOT=   zeros(N,1);
w=      zeros(N,Y);
C=      zeros(N,Y);
E=      zeros(N,Y);
E_TOT=  zeros(N,1);
Estar=  zeros(N,1);
HH=     zeros(N,Y);
H_TOT=  zeros(N,1);
S=      zeros(N,Y);
S_t=    zeros(1,Y);
S_r=    zeros(1,Y);
S_v=    zeros(1,Y);
S_el=   zeros(1,Y);
S_TOT=  zeros(N,1);
Sstar=  zeros(N,1);
cv_tot= zeros(N,1);
cp_tot= zeros(N,1);
Z=      zeros(N,1);
eel=    zeros(1,Y);
ecomp=  zeros(N,1);
scomp=  zeros(N,1);

% Zero energy state corrections
e0(1) = 0;
e0(2) = .5*R(2)*TH_D(1);
e0(3) = 0;
e0(4) = .5*R(4)*TH_D(3);
e0(5) = R(5)*(-TH_D(5) + .5*(TH_D(1)+TH_D(3)));
e0(6) = 0;
e0(7) = R(7)*(TH_D(3)/2 + TH_I(4));
e0(8) = R(8)*(TH_D(1)/2 + TH_I(2));
%% Solution
for o = 1:1
    Pin = Pspan(o);
    P = Pin*101325; %PA
    PARTPi = zeros(1,Y)+P/Y;
    rho0 = P/(R0*300);
    % CALCULATE SPECIES RELATIONS
    for i = 1:1
        T(i) = Tmin + z*(i-1);
        
        % PARTITION FUNCTIONS
        Q_el(1) = 1;
        Q_el(2) = 4 + 10*exp(-27700/T(i)) + 6*exp(-41500/T(i));
        Q_el(3) = 3 + 2*exp(-11390/T(i)) + 2*exp(-18990/T(i));
        Q_el(4) = 5 + 3*exp(-228/T(i)) + 1*exp(-326/T(i)) ...
            + 5*exp(-22800/T(i)) + 1*exp(-48600/T(i));
        Q_el(5) = 2 + 2*exp(-174/T(i));
        Q_el(6) = 2;
        Q_el(7) = 4 + 10*exp(-38600/T(i)) + 6*exp(-58200/T(i));
        Q_el(8) = 1 + 3*exp(-70.6/T(i)) + 5*exp(-188.9/T(i)) + ...
            5*exp(-22000/T(i)) + 1*exp(-47000/T(i)) + 5*exp(-67900/T(i));
        for j =1:Y
            Q_t(j) = (2*pi*M(j)*K*T(i)/(H^2))^(3/2);
            if j == 1 || j == 3 || j == 5
                Q_r(j) = T(i)/(SIGMA(j)*TH_R(j));
                Q_v(j) = 1/(1-exp(-TH_V(j)/T(i)));
            else
                Q_r(j) = 1;
                Q_v(j) = 1;
            end
            Q(i,j) = Q_t(j)*Q_r(j)*Q_v(j)*Q_el(j);
        end
        
        KP(i,1) = 1/(K*T(i)) * exp(TH_D(1)/T(i)) * Q(i,1)/Q(i,2)^2;
        KP(i,2) = 1/(K*T(i)) * exp(TH_D(3)/T(i)) * Q(i,3)/Q(i,4)^2;
        KP(i,3) = 1/(K*T(i)) * exp(TH_D(5)/T(i)) * Q(i,5)/(Q(i,2)*Q(i,4));
        KP(i,4) = 1/(K*T(i)) * exp(TH_I(4)/T(i)) * Q(i,4)/(Q(i,6)*Q(i,7));
        KP(i,5) = 1/(K*T(i)) * exp(TH_I(2)/T(i)) * Q(i,2)/(Q(i,6)*Q(i,8));
        
        kpi = KP(i,:);
        options=optimset('Display','off');
        PARTP(i,:) = fsolve(@(PARTPi)eqns(PARTPi,kpi,P),PARTPi,options);
        PARTPi(:) = PARTP(i,:);
        RHOT(i)=0;
        P_TOT(i)=0;
        for j=1:Y
            P_TOT(i) = PARTP(j)+P_TOT(i);
            RHO(i,j) = PARTP(i,j)/(R(j)*T(i)); %kg/m3 DENSITY
            RHOT(i) = RHOT(i) + RHO(i,j);
            X(i,j)=PARTP(i,j)/(RGAS*T(i)); %kmol/m3 MOLE CONCENTRATION
        end
        
        %Calculate mass fraction and mole fraction
        temp=0;
        for j=1:Y
            rhostar(i,j) = RHO(i,j)/rho0;
            C(i,j)=RHO(i,j)/RHOT(i); % MASS FRACTION
            temp=temp+C(i,j)/M(j); % USED TO FIND MASS-AVERAGED MOLECULAR WEIGHT
        end
        mMT=1/temp;
        for j=1:Y
            w(i,j) = C(i,j)*mMT/M(j); % MOLE FRACTION
        end
    end
    
    % FIND COMPRESSIBILITY FACTOR
    for i =1:N
        Z(i)=P/(RHOT(i)*R0*T(i));
    end
    
    % FIND SPECIFIC INTERNAL ENERGY & DIMENSIONLESS INTERNAL ENERGY
    % COMMON REFERENCE STATE OF ZERO ENERGY IS SET AT THE COMPLETE COMBINED
    % DIATOMIC STATE
    for i = 1:N
        
        for j = 1:Y
            eel(j) = R(j)*th_el(2,j)*(g(2,j)/g(1,j))* ...
                exp(-th_el(2,j)/T(i))/(1+(g(2,j)/g(1,j))*exp(-th_el(2,j)/T(i)));
        end
        
        E(i,1) = 5/2*R(1)*T(i) + R(1)*TH_V(1)/(exp(TH_V(1)/T(i))-1) + eel(1) + e0(1); %J/KG
        E(i,3) = 5/2*R(3)*T(i) + R(3)*TH_V(3)/(exp(TH_V(3)/T(i))-1) + eel(2) + e0(3); %J/KG
        E(i,5) = 5/2*R(5)*T(i) + R(5)*TH_V(5)/(exp(TH_V(5)/T(i))-1) + eel(3) + e0(5); %J/KG
        E(i,2) = 3/2*R(2)*T(i) + eel(4) + e0(2); %J/KG
        E(i,4) = 3/2*R(4)*T(i) + eel(5) + e0(4); %J/KG
        E(i,6) = 3/2*R(6)*T(i) + eel(6) + e0(6); %J/KG
        E(i,7) = 3/2*R(7)*T(i) + eel(7) + e0(7); %J/KG
        E(i,8) = 3/2*R(8)*T(i) + eel(8) + e0(8); %J/KG
        
        E_TOT(i) = 0; %J/KG
        RT(i) = 0;
        for j=1:Y
            e_temp =C(i,j)*E(i,j);
            r_temp =C(i,j)*R(j);
            E_TOT(i)=E_TOT(i)+e_temp;
            RT(i) = RT(i) + r_temp;
        end
        
        % Dimensionless specific internal energy
        Estar(i) = Z(i)*E_TOT(i)/(RT(i)*T(i));
        
        %Comparison of electronic and translational energy modes
        ecomp(i) = 0; %J/KG
        for j=1:Y
            ecomp_temp =C(i,j)*eel(j)/(1.5*R(j)*T(i));
            ecomp(i)=ecomp(i)+ecomp_temp;
        end
    end
    
    % FIND SPECIFIC ENTHALPY
    for i=1:N
        H_TOT(i)=0; %J/KG
        for j=1:Y
            HH(i,j) = E(i,j) + R(j)*T(i); % J/KG
            h_temp =C(i,j)*HH(i,j);
            H_TOT(i)=H_TOT(i)+h_temp;
        end
        hstar(i) = Z(i)*H_TOT(i)/(RT(i)*T(i));
    end
    
    % FIND SPECIFIC ENTROPY & DIMENSIONLESS ENTROPY
    for i = 1:N
        for j = 1:Y
            S_el(j) = R(j)*(log(g(1,j)) + log(1+(g(2,j)/g(1,j))*exp(-th_el(2,j)/T(i))) + ...
                (g(2,j)/g(1,j))*(th_el(2,j)/T(i))*exp(-th_el(2,j)/T(i))/(1+(g(2,j)/g(1,j))*exp(-th_el(2,j)/T(i))));
        end
        for j =1:Y
            S_t(j) = 5/2*R(j)*log(T(i)) - R(j)*log(P) + R(j)*(log((2*pi*M(j)/(H^2))^(3/2)*K^(5/2))+5/2);
            if j == 1 || j == 3 || j == 5
                S_v(j) = R(j)*(-log(1-exp(-TH_V(j)/T(i))) + TH_V(j)/T(i)/(exp(TH_V(j)/T(i))-1));
                S_r(j) = R(j)*(log(T(i)/(SIGMA(j)*TH_R(j)))+1);
            else
                S_r(j) = 0;
                S_v(j) = 0;
            end
            S(i,j) = S_t(j)+S_r(j)+S_v(j)+S_el(j);
        end
        
        S_TOT(i) = 0; %J/KG-K
        for j=1:Y
            s_temp =C(i,j)*S(i,j);
            S_TOT(i)=S_TOT(i)+s_temp;
        end
        % Dimensionless specific entropy
        Sstar(i) = Z(i)*S_TOT(i)/(RT(i));
        
        scomp(i) = 0; %J/KG
        for j=1:Y
            scomp_temp =C(i,j)*S_el(j)/S_t(j);
            scomp(i)=scomp(i)+scomp_temp;
        end
    end
    
    % Find specific heat
    for i =1:N
        for j = 1:Y
            cv_t(j) = (3/2)*R(j);
            %cv_el(j) = 0;
            cv_el(j) = R(j)*(th_el(2,j)/T(i))^2*(g(2,j)/g(1,j))*exp(-th_el(2,j)/T(i))/(1+(g(2,j)/g(1,j))*exp(-th_el(2,j)/T(i)))^2;
            if j == 1 || j == 3 || j == 5
                cv_v(j) = R(j)*((TH_V(j)/(2*T(i)))/(sinh(TH_V(j)/(2*T(i)))))^2;
                cv_r(j) = R(j);
            else
                cv_r(j) = 0;
                cv_v(j) = 0;
            end
            cv(i,j) = cv_t(j)+cv_r(j)+cv_v(j)+cv_el(j);
        end
        cv_tot(i)=0;
        for j=1:Y
            cv_temp =C(i,j)*cv(i,j);
            cv_tot(i)=cv_tot(i)+cv_temp;
        end
        cp_tot(i) = cv_tot(i) + RT(i);
        % Dimensionless specific entropy
        cstar(i) = Z(i)*cv_tot(i)/(RT(i));
        gamma(i) = cp_tot(i)/cv_tot(i);
        a_ideal = sqrt(1.4*R0*T(i));
        a(i) = sqrt(gamma(i)*RT(i)*T(i));
        astar(i) = a(i)/a_ideal;
        
    end
    
    
    cd('DATA')
    fOut=sprintf('Estar_%03d.mat',Pin);
    save(fOut,'Estar')
    fOut=sprintf('Sstar_%03d.mat',Pin);
    save(fOut,'Sstar')
    fOut=sprintf('hstar_%03d.mat',Pin);
    save(fOut,'hstar')
    fOut=sprintf('cstar_%03d.mat',Pin);
    save(fOut,'cstar')
    fOut=sprintf('gamma_%03d.mat',Pin);
    save(fOut,'gamma')
    fOut=sprintf('Z_%03d.mat',Pin);
    save(fOut,'Z')
    fOut=sprintf('w_%03d.mat',Pin);
    save(fOut,'w')
    fOut=sprintf('rhostar_%03d.mat',Pin);
    save(fOut,'rhostar')
    fOut=sprintf('T_%03d.mat',Pin);
    save(fOut,'T')
    fOut=sprintf('ecomp_%03d.mat',Pin);
    save(fOut,'ecomp')
    fOut=sprintf('astar_%03d.mat',Pin);
    save(fOut,'astar')
    cd('../')
    
    % TEMPERATURE AND MOLE FRACTION
   % myColors={'r','k','b','g','m','c','y','r','k'};
    %figure(o)
    %hold on
    for j=1:Y
        %plot(T,w(:,j),'Color',myColors{j});
        %plot(T,C(:,j),'Color',myColors{j});
        %plot(T,RHO(:,j),'Color',myColors{j});
    end
    %legend('N_2','N','O_2','O','NO','e^-','O^+','N^+')
    %ylabel('w - Mole Fraction')
    %ylabel('C - Mass fraction')
    %ylabel('\rho^* (\rho/\rho_0) - density [kg/m^3]')
    %xlabel('T - Temperature [K]')
    %xlim([2000 15000])
    %ylim([0 1])
    
end
