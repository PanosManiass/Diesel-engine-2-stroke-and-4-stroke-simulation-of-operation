

%------------------- Panagiotis Manias ------------------------------------

clc
clear all
close all
    %-------------------------Default Input Data---------------------

%----------Engine Inputs------
V=6311                      %Engine Displacement Volume(Liter)
N=8;                      %Number of Cylinders
n=100;                   %Engine Speed(rpm)
r_c=32;                   %Compression Ratio
X=4;                   %Stroke-to-Bore Ratio S/B
eta_m=0.89;               %Mechanical Efficiency
EGR=5;                    %Internal EGR Percentage
P_Turbo=110;               %Turbo Boost
P_Exhaust=120;            %Exhaust Pressure

%---------Fuel Inputs/Combustion Efficiency----

x=1;                     %number of Carbon atoms
y=4;                     %number of Hydrogen atoms
AF=17.2;                    %Air-Fuel Mass Ratio
eta_conv=0.52;               %Combustion Efficiency
MR_CH4=16.4                   %gr/mol
MR_CO2=44.0                    %gr/mol

%-----------Properties of Working Fluids----------------
%a constant average value between two exteremes ( SOC & EOCombustion)
K=1.33;                   %Specific Heat Ratio of Air
R=0.287;                  %Particular Gas Constant of Air                (kJ/kg.K)
R_NG=0.5183                   %Particular Gas Constant of NG                (kJ/kg.K)
C_pNG=4.595                    %Specific Heat for Constant Pressure of NG @ 1000K    (kJ/kg)
C_v=0.821;                %Specific Heat for Constant Volume of Air      (kJ/kg)
C_p=R+C_v;                %Specific Heat for Constant Pressure of Air    (kJ/kg)

%---------Atmospheric Inputs------------
P_0=101;                  %Atmospheric Pressure                                     (kPa)
T_0=45;                   %Intake Temperature                (C)
row_0=P_0/(R*(293.15));   %Density of Air at Reference Condition                    (kg/cubic meter)     

%--------------------------Calculation---------------------------


%--------Fuel Heating Value & Combustion Efficiency----
FHV=50000;               %Fuel Heating Value (kJ/kg) [Dulong’s formula]b=x/y;
AF_s=17.2;             %Stoichiometric Air/Fuel ratio (Heywood (5-3))
lambda=AF/AF_s;
Th_Air=lambda*100;                                 %Theoritical Air Percentage
EAP=Th_Air-100;                                    %Exess Air Percentage

if lambda >= 1                                     %Schmidt (1996)
eta_ch=1;
else
    eta_ch=(1.3773*lambda)-0.3773;
end

eta_c=eta_conv*eta_ch;                             %Pischinger (1989)
%-----------------------------------------------------------------

%for one cylinder
V_d=(V/N)/1000;                 %Displacement Volume (Cubic Meter)
V_c=V_d/(r_c-1);                %Clearance Volume    (Cubic Meter)
B=((4*V_d/pi*X)^(1/3))*100;     %Bore                (cm)
S=B*X;                          %Stroke              (cm)
A_p=(pi*(B^2))/4 ;              %Piston Area        (square cm)

%----State 1 (Start of Compression)-----
T_1=(T_0)+273;                  %Temprature of Gas Mixture in Cylinder at State 1 (Deg K)
P_1=P_0 + P_Turbo;                        %Pressure of Gas Mixture in Cylinder at State 1   (kPa)
V_1=V_d+V_c;                    %Volume of Gas Mixture in Cylinder at State 1     (Cubic Meter)

%----Mass of Gases Calculation------
M_m=(P_1*V_1)/(R*T_1);   %Mass of Gas Mixture in Cylinder (kg)

%the mass within the cylinder will then remain the same for the entire
%cycle

M_a=((100-EGR)/100)*M_m;    %Mass of Air in Cylinder                         (kg)
M_f=(1/((1.95)*(AF)))*M_a;     %Mass of Fuel in Cylinder              Excess air principle due to direct injection, lamda=2          (kg)
M_ex=(EGR/100)*(M_f+M_a);                     %Mass of Exhaust Residual in Cylinder at SOC     (kg)

%------------State 2 (Start of Combustion)-----

T_2=(T_1)*((r_c)^(K-1));                %Temprature of Gas Mixture in Cylinder at State 2 (Deg K)
P_2=(P_1)*((r_c)^(K));                    %Pressure of Gas Mixture in Cylinder at State 2   (kPa)
V_2=V_c;                            %Volume of Gas Mixture in Cylinder at State 2     (Cubic Meter)


%------------State 3 (End of Combustion)-----

T_3=((FHV*eta_c)/((AF+1)*(C_p+0.047)))+T_2;     %Temprature of Gas Mixture in Cylinder at State 3 (Deg K)
P_3=P_2;                                %Pressure of Gas Mixture in Cylinder at State 3   (kPa)
V_3=(T_3/T_2)*V_2;                      %Volume of Gas Mixture in Cylinder at State 3     (Cubic Meter)

%----Cutoff Ratio (Change in Volume During Combustion)----
Beta=V_3/V_2;

%------------State 4 (End of Expansion)-----
V_4=V_1;                                %Volume of Gas Mixture in Cylinder at State 4     (Cubic Meter)
T_4=T_3*((V_3/V_4)^(K-1));              %Temprature of Gas Mixture in Cylinder at State 4 (Deg K)
P_4=P_3*((V_3/V_4)^(K));                  %Pressure of Gas Mixture in Cylinder at State 4   (kPa)


%------------------------------------------------------------------------------------------------------

%---------------------------------Processes------------------------------
%--- 1-2 Process (Compression)
Q_12=0;
W_12=M_m*C_v*(T_1-T_2);     %(kJ)

%--- 2-3 Process (Constant-Pressure Combustion)
Q_in=(M_f*C_pNG+M_m*C_p)*(T_3-T_2);     %(kJ)
W_23=P_3*(V_3-V_2);         %(kJ)

%--- 3-4 Process (Expansion)
Q_34=0;
W_34=M_m*C_v*(T_3-T_4);     %(kJ)

%--- 4-5 Process  (Exhaust)
Q_41=M_m*C_v*(T_1-T_4);     %(kJ)   
W_41=0;


%------------------------------ Overal Engine Performance----------

U_p=S*n*(1/6000);                         %Mean Piston Speed                                    (m/s)
W_i=W_12+W_23+W_34;                         %Net Indicated Work for One Cylinder During One Cycle (kJ)
eta_t=100*W_i/Q_in;                             %Thermal Efficiency
imep=W_i/(V_d);                             %Indicated Mean Effective Pressure                    (kPa)
bmep=imep*eta_m;                            %Brake Mean Effective Pressure                        (kPa)
W_idat=W_i*N*(n/(60));                    %Indicated Power at {n} RPM                           (KW)
W_b=eta_m*W_i;                              %Net Brake Work for One Cylinder During One Cycle     (kJ)
W_bdat=W_b*N*(n/(60));                    %Brake Power at {n} RPM                               (KW)
Taw=(W_bdat/(2*pi*(n/60)))*1000;            %Torque                                               (N-m)
W_fdat=W_idat-W_bdat;                       %Friction Power Lost                                  (KW)
BSP=W_bdat/A_p;                             %Brake Specific Power                                 (KW/square cm)
OPD=W_bdat/V;                               %Output per Displacement                              (KW/Litr)
M_fdat=M_f*N*1000*n*60;                  %Fuel That is Injected in All Cylinders During One hour(g/cycle)
bsfc=(M_fdat/W_bdat);               %Brake Specific Fuel Consumption                      (gr/KW-hr)
eta_v=M_a/(row_0*V_d);                      %Volumetric Efficiency
CO2_Emission=bsfc*MR_CO2/MR_CH4


%--------------------Display Outputs--------------------------

disp('                                                                   ');
disp('                                                                   ');

disp('---------------------------Fuel & Combustion Properties-------------------------- ');

F=[x y];
disp('     C     H')
disp(F)
A = ['The Actual Air/Fuel Ratio is: ',num2str(AF)];
    disp(A)
    
Z = ['The Stoichiometric Air/Fuel Ratio is: ',num2str(AF_s)];
disp(Z)

T=['The Percent of Theoritical Air is: ',num2str(Th_Air)];
disp(T)

E=['The Percent of Excess Air is: ',num2str(EAP)];
disp(E)

Eff=['Total Combustion Efficiency is: ',num2str(eta_c)];
disp(Eff)

HV=['Fuel Heating Value is: ',num2str(FHV), ' kJ/kg'];
disp(HV)

disp('                                                                   ');
disp('                                                                   ');
disp('---------------------------Cycle Operating Point Properties-------------------------- ');
disp('                                                                   ');
disp('                                                                   ');
aa=['The Maximum Temprature of The Cycle is: ',num2str(T_3-273.15), ' Deg C'];
    disp(aa)

ab=['The Maximum Pressure of The Cycle is : ',num2str(P_3), ' kPa'];
    disp(ab)
disp('                                                                   ');
disp('                                                                   ');

disp('---------------------------Overal Engine Performances-------------------------- ');
disp('                                                                   ');
disp('                                                                   ');

a15=['The Mean Piston Speed is: ',num2str(U_p),' m/s'];
    disp(a15)


a1=['The Thermal Efficiency is: ',num2str(eta_t), ' %'];
    disp(a1)

a2=['The Indicated Mean Effective Pressure imep is : ',num2str(imep), ' kPa'];
    disp(a2)

a3 = ['The Brake Mean Effective Pressure is : ',num2str(bmep), ' kPa'];
    disp(a3)
    
a5 = ['The Indicated Power at ',num2str(n),' RPM is : ',num2str(W_idat), ' KW'];
    disp(a5)

a6 = ['The Brake Power at ',num2str(n),' RPM is : ',num2str(W_bdat), ' KW'];
    disp(a6)


a8=['The Friction Power Lost is: ',num2str(W_fdat), ' KW'];
    disp(a8)

    
a7 = ['The Torque at ',num2str(n),' RPM is : ',num2str(Taw), ' N.m'];
    disp(a7)    
 
    
a9=['The Brake Specific Fuel Consumption bsfc is : ',num2str(bsfc), ' gr/KW-hr'];
    disp(a9)

a12 = ['The Volumetric Efficiency is : ',num2str(eta_v)];
    disp(a12)
        
a10 = ['The Output per Displacement OPD is : ',num2str(OPD), ' KW/Litr'];
    disp(a10)

a11 = ['The Brake Specific Power BSP is : ',num2str(BSP), ' KW/cm^2'];
    disp(a11)

a14 = ['The CutOff Ratio is : ',num2str(Beta)];
    disp(a14)

a15 = ['The CO2 emissions per kWh are : ',num2str(CO2_Emission), ' gr.'];
    disp(a15)


%-------------------%P-V Diagramm----------------------

%---Plot 1-2 Process (Compression)
V_12=linspace(V_1,V_2,180);
P_12=((V_1./V_12).^K)*P_1;
plot(V_12,P_12);
title('P-V Diagram','FontSize', 20 , 'FontName', 'Helvetica','FontWeight','bold');xlabel('Volume (m^3)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');ylabel('Pressure (kPa)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');
hold on

%---Plot 2-3 Process (Combustion)
V_23=[V_2 V_3];
P_23=[P_2 P_3];
plot(V_23,P_23);

%---Plot 3-4 Process (Expansion)
V_34=linspace(V_3,V_4,180);
P_34=((V_3./V_34).^K)*P_3;
plot(V_34,P_34);

%---Plot 4-1 Process  (Exhaust)
V_41=[V_4 V_1];
P_41=[P_4 P_1];
plot(V_41,P_41);

hold off

text(V_1,P_1, '1' ,  'FontSize' ,12)
text(V_2,P_2, '2' ,  'FontSize' ,12)
text(V_3,P_3, '3' ,  'FontSize' ,12)
text(V_4,P_4, '4' ,  'FontSize' ,12)


figure;

%---------------- T-S Diagram-----

%Refference State:   ( Table A/12SI Van Wylen Ed.4 )
Pref=P_0;            % 0.1 (MPa)
Tref=298.15;         % 25 Deg C
Sref=6.86305;        % (kJ/kgK)

%---Plot 1-2 Process (Compression)

S_1=Sref+(C_p*log((T_0+273.15)/Tref));      %( Eq.10-31 Van Wylen Ed.4 )
S_2=S_1;
S_12(1:1, 1:180)=S_2;
T_12=T_1*((V_1./V_12).^(K-1));
plot(S_12,T_12);
title('T-S Diagram','FontSize', 20 , 'FontName', 'Helvetica','FontWeight','bold');xlabel('Entropy (kJ/kg)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');ylabel('Temprature (K)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');
hold on

%---Plot 2-3 Process (Combustion)
T_23=linspace(T_2,T_3,32);
S_23=S_2+(C_p*log(T_23/T_2));           %( Eq.10-33 Van Wylen Ed.4 )
S_3=S_2+(C_p*log(T_3/T_2));
plot(S_23,T_23);

%---Plot 3-4 Process (Expansion)
S_4=S_3;
S_34(1:1, 1:180)=S_4;
T_34=T_3*((V_3./V_34).^(K-1));
plot(S_34,T_34);

%---Plot 4-1 Process  (Exhaust)
T_41=linspace(T_4,T_1,180);
S_41=S_4+(C_v*log(T_41/T_4));           %( Eq.10-33 Van Wylen Ed.4 )
plot(S_41,T_41);


hold off

text(S_1,T_1, '1' ,  'FontSize' ,12)
text(S_2,T_2, '2' ,  'FontSize' ,12)
text(S_3,T_3, '3' ,  'FontSize' ,12)
text(S_4,T_4, '4' ,  'FontSize' ,12)