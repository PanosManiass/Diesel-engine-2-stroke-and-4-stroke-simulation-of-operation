%------------------- Panagiotis Manias ------------------------------------
clc
clear all
close all
    %-------------------------Default Input Data---------------------
%----------Engine Inputs------
V=461.7;                    %Engine Displacement Volume(Liters)
N=10;                      %Number of Cylinders
n=750;                   %Engine Speed(rpm)
r_c=9.5;                  %Compression Ratio
X=1.257;                  %Stroke-to-Bore Ratio S/B
eta_m=0.90;               %Mechanical Efficiency
EGR=15;                 %Internal EGR Percentage
P_Turbo=100;               %Turbo Boost
P_Exhaust=105;            %Exhaust Pressure
Diesel_Percentage = 0.4
Ammonia_Percentage = 0.6
Hydrogen_Percentage = 0.0

%---------Fuel Inputs/Combustion Efficiency----

x=12;                     %number of Carbon atoms
y=23;                     %number of Hydrogen atoms
AF_Diesel=20;                    %Air-Fuel Mass Ratio of Diesel
AF_NH3=6.0;                    %Air-Fuel Mass Ratio of Ammonia
AF_H2=34.0;                    %Air-Fuel Mass Ratio of Hydrogen
eta_conv=0.5;               %Combustion Efficiency
MR_Diesel = x*12+y*1            %gr/mol
MR_H2=2                         %gr/mol
MR_NH3=17.031                   %gr/mol
MR_CO2=44.0                     %gr/mol


%Properties of Working Fluids---------------- 
%a constant average value between two exteremes ( SOC & EOCombustion)
K=1.33;                   %Specific Heat Ratio of Air
R=0.287;                  %Particular Gas Constant of Air                (kJ/kg.K)
R_NG=0.5183               %Particular Gas Constant of NG                  (kJ/kg.K)
C_vNG=1.7354               %Specific Heat for Constant Pressure of NG @ 1000K    (kJ/kg)
C_v=0.821;                %Specific Heat for Constant Volume of Air      (kJ/kg)
C_p=R+C_v;                %Specific Heat for Constant Pressure of Air    (kJ/kg)

%---------Atmospheric Inputs------------
P_0=101;                  %Atmospheric Pressure                                  (kPa)
T_0=40;                   %Combustion Chamber Temperature at Start of Compression
row_0=P_0/(R*(293.15));   %Density of Air at Reference Condition                 (kg/cubic meter) 

%--------------------------Calculation---------------------------

%--------Fuel Heating Value & Combustion Efficiency----
FHV_Diesel=42800;                 %Fuel Heating Value (kJ/kg)
FHV_Hydrogen=120000;               %Fuel Heating Value (kJ/kg)
FHV_Ammonia=18500;   
AF_s=17.2;             %Stoichiometric Air/Fuel ratio (Heywood (5-3))
lambda=AF/AF_s;
Th_Air=lambda*100;                                 %Theoritical Air Percentage
EAP=Th_Air-100;                                    %Exess Air Percentage

if lambda >= 1                                    
eta_ch=1;
else
    eta_ch=(1.3773*lambda)-0.3773;                 %Schmidt (1996)
end

eta_c=eta_conv*eta_ch;                             %Pischinger (1989)
%-----------------------------------------------------------------


%for one cylinder
V_d=(V/N)/1000;                 %Displacement Volume (Cubic Meter)
V_c=V_d/(r_c-1);                %Clearance Volume    (Cubic Meter)
B=((4*V_d/pi*X)^(1/3))*100;     %Bore                (cm)
S=B*X;                          %Stroke              (cm)
A_p=(pi*(B^2))/4 ;               %Piston Area         (square cm)
%----State 1 (Start of Compression)-----
T_1=(T_0)+273;                  %Temprature of Gas Mixture in Cylinder at State 1 (Deg K)
P_1=P_Turbo+P_0;                        %Pressure of Gas Mixture in Cylinder at State 1   (kPa)
V_1=V_d+V_c;                    %Volume of Gas Mixture in Cylinder at State 1     (Cubic Meter)

%----Mass of Gases Calculation------
M_m=(P_1*V_1)/(R*T_1);   %Mass of Gas Mixture in Cylinder (kg)

%the mass within the cylinder will then remain the same for the entire
%cycle

M_a=((100-EGR)/100)*M_m;    %Mass of Air in Cylinder                         (kg)
M_f=((Diesel_Percentage*1/((1.3)*(AF_Diesel)))+(Ammonia_Percentage*1/((1.95)*(AF_NH3)))+(Hydrogen_Percentage*1/((1.95)*(AF_H2))))*M_m;     %Mass of Fuel in Cylinder                        (kg)
M_ex=(EGR/100)*(M_f+M_a);                     %Mass of Exhaust Residual in Cylinder at SOC     (kg)
   
%------------State 2 (Start of Combustion)-----

T_2=(T_1)*((r_c)^(K-1));                %Temprature of Gas Mixture in Cylinder at State 2 (Deg K)
P_2=(P_1)*((r_c)^K);                    %Pressure of Gas Mixture in Cylinder at State 2   (kPa)
V_2=V_1/r_c;                            %Volume of Gas Mixture in Cylinder at State 2     (Cubic Meter)


%------------State 3 (End of Combustion)-----

T_3=((FHV*eta_c)/((AF)*(C_v)))+T_2;     %Temprature of Gas Mixture in Cylinder at State 3 (Deg K)
P_3=(P_2/T_2)*T_3;                      %Pressure of Gas Mixture in Cylinder at State 3   (kPa)
V_3=V_2;                                %Volume of Gas Mixture in Cylinder at State 3     (Cubic Meter)

%----Pressure Ratio (Change in Pressure During Combustion)----
Alpha=P_3/P_2;                          %Pressure Increase During Contant-Volume Combustion 

%------------State 4 (End of Expression)-----
T_4=T_3*((1/r_c)^(K-1));                %Temprature of Gas Mixture in Cylinder at State 4 (Deg K)
P_4=P_3*((1/r_c)^K);                    %Pressure of Gas Mixture in Cylinder at State 4   (kPa)
V_4=(M_m*R*T_4)/P_4;                    %Volume of Gas Mixture in Cylinder at State 4     (Cubic Meter)

%------------State 5 
T_5=T_1;                                %Temprature of Gas Mixture in Cylinder at State 5 (Deg K)
P_5=P_1;                                %Pressure of Gas Mixture in Cylinder at State 5   (kPa)
V_5=V_1;                                %Volume of Gas Mixture in Cylinder at State 5     (Cubic Meter)

%------------State 6 
T_6=T_0;                                %Temprature of Gas Mixture in Cylinder at State 6 (Deg K)
P_6=P_0+20;                                %Pressure of Gas Mixture in Cylinder at State 6   (kPa)
V_6=V_2;                                %Volume of Gas Mixture in Cylinder at State 6     (Cubic Meter)

%------------------------------------------------------------------------------------------------------

%---------------------------------Processes------------------------------
%--- 1-2 Process (Compression)
Q_12=0;
W_12=M_m*C_v*(T_1-T_2);     %(kJ)

%--- 2-3 Process (Combustion)
Q_23=M_m*FHV*eta_c;         %(kJ)
W_23=0;

%--- 3-4 Process (Expansion)
Q_34=0;
W_34=M_m*C_v*(T_3-T_4);     %(kJ)

%--- 4-5 Process  (Exhaust)
Q_45=M_ex*C_v*(T_5-T_4);     %(kJ)   
W_45=0;

%------------------------------ Overal Engine Performance----------

U_p=2*S*n*(1/6000);                         %Mean Piston Speed                                    (m/s)
W_i=W_12+W_34;                              %Net Indicated Work for One Cylinder During One Cycle (kJ)
Q_in=M_f*FHV;                         %Heat Added for One Cylinder During One Cycle         (kJ)
eta_t=100*(W_i)/Q_in;                             %Thermal Efficiency
imep=W_i/(V_d);                             %Indicated Mean Effective Pressure                    (kPa)
bmep=imep*eta_m;                            %Brake Mean Effective Pressure                        (kPa)
W_idat=W_i*N*(n/(60*2));                    %Indicated Power at {n} RPM                           (KW)
W_b=eta_m*W_i;                              %Net Brake Work for One Cylinder During One Cycle     (kJ)
W_bdat=W_b*N*(n/(60*2));                    %Brake Power at {n} RPM                               (KW)
Taw=(W_bdat/(2*pi*(n/60)))*1000;            %Torque                                               (N-m)
W_fdat=W_idat-W_bdat;                       %Friction Power Lost                                  (KW)
BSP=W_bdat/A_p;                             %Brake Specific Power                                 (KW/square cm)
OPD=W_bdat/V;                               %Output per Displacement                              (KW/Litr)
M_fdat=M_f*(n/60)*(1/2)*N;                  %Rate of Fuel That Inject in All Cylinder During One Cycle (kg/sec)
bsfc=(M_fdat/W_bdat)*3600000;               %Brake Specific Fuel Consumption                      (gr/KW-hr)
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

    
a7 = ['The Torque at ',num2str(n),' RPM is : ',num2str(Taw), ' N-m'];
    disp(a7)    
 
    
a9=['The Brake Specific Fuel Consumption bsfc is : ',num2str(bsfc), ' gr/KW-hr'];
    disp(a9)

a12 = ['The Volumetric Efficiency is : ',num2str(eta_v)];
    disp(a12)
        
a10 = ['The Output per Displacement OPD is : ',num2str(OPD), ' KW/Litr'];
    disp(a10)

a11 = ['The Brake Specific Power BSP is : ',num2str(BSP), ' KW/cm^2'];
    disp(a11)

a13 = ['The Pressure Increase Ratio During Contant-Volume Combustion is : ',num2str(Alpha)];
    disp(a13)

a15 = ['The CO2 emissions per kWh are : ',num2str(CO2_Emission), ' gr.'];
    disp(a15)
 
%-------------------------Plot Statements-------------------------

%---------------- P-V Diagram-----

%---Plot 1-2 Process (Compression)

V_12=linspace(V_1,V_2,100);
P_12=((V_1./V_12).^K)*P_1;
plot(V_12,P_12);
title('P-V Diagram','FontSize', 20 , 'FontName', 'Helvetica','FontWeight','bold');xlabel('Volume (m^3)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');ylabel('Pressure (kPa)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');
hold on

%---Plot 2-3 Process (Combustion)
V_23=[V_2 V_3];
P_23=[P_2 P_3];
plot(V_23,P_23);

%---Plot 3-4 Process (Expansion)
V_34=linspace(V_3,V_4,100);
P_34=((V_3./V_34).^K)*P_3;
plot(V_34,P_34);

%---Plot 4-1 Process  (Exhaust)
V_41=[V_4 V_1];
P_41=[P_4 P_1];
plot(V_41,P_41);

%---Plot 5-6 Process (Pumping-Expulsion)
V_56=[V_5 V_6];
P_56=[P_5 P_6];
plot(V_56,P_56);


hold off

text(V_1,P_1, '1,5' ,  'FontSize' ,12)
text(V_2,P_2, '2' ,  'FontSize' ,12)
text(V_3,P_3, '3' ,  'FontSize' ,12)
text(V_4,P_4, '4' ,  'FontSize' ,12)
text(V_6,P_6, '6' ,  'FontSize' ,12)

figure;

%---------------- T-S Diagram-----

%Refference State:   ( Table A/12SI Van Wylen Ed.4 )
Pref=P_0;            % 0.1 (MPa)
Tref=298.15;         % 25 Deg C
Sref=6.86305;        % (kJ/kgK)

%---Plot 1-2 Process (Compression)

S_1=Sref+(C_p*log((T_0+273.15)/Tref));      %( Eq.10-31 Van Wylen Ed.4 )
S_2=S_1;
S_12=[S_1 S_2];
T_12=[T_1 T_2];
plot(S_12,T_12);
title('T-S Diagram','FontSize', 20 , 'FontName', 'Helvetica','FontWeight','bold');xlabel('Entropy (kJ/kg)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');ylabel('Temprature (K)','FontSize', 14 , 'FontName', 'Helvetica','FontWeight','bold');
hold on

%---Plot 2-3 Process (Combustion)
T_23=linspace(T_2,T_3,100);
S_23=S_2+(C_v*log(T_23/T_2));           
S_3=S_2+(C_v*log(T_3/T_2));
plot(S_23,T_23);

%---Plot 3-4 Process (Expansion)
S_4=S_3;
S_34=[S_3 S_4];
T_34=[T_3 T_4];
plot(S_34,T_34);

%---Plot 4-1 Process  (Exhaust)
T_41=linspace(T_4,T_1,100);
S_41=S_4+(C_v*log(T_41/T_4));           
plot(S_41,T_41);


hold off

text(S_1,T_1, '1' ,  'FontSize' ,12)
text(S_2,T_2, '2' ,  'FontSize' ,12)
text(S_3,T_3, '3' ,  'FontSize' ,12)
text(S_4,T_4, '4' ,  'FontSize' ,12)