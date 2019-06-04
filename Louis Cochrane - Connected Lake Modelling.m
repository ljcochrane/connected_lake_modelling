% 
%    ____                            _           _   _          _             
%   / ___|___  _ __  _ __   ___  ___| |_ ___  __| | | |    __ _| | _____  ___ 
%  | |   / _ \| '_ \| '_ \ / _ \/ __| __/ _ \/ _` | | |   / _` | |/ / _ \/ __|
%  | |__| (_) | | | | | | |  __/ (__| ||  __/ (_| | | |__| (_| |   <  __/\__ \
%   \____\___/|_| |_|_| |_|\___|\___|\__\___|\__,_| |_____\__,_|_|\_\___||___/
%  |  \/  | ___   __| | ___| | (_)_ __   __ _                                 
%  | |\/| |/ _ \ / _` |/ _ \ | | | '_ \ / _` |                                
%  | |  | | (_) | (_| |  __/ | | | | | | (_| |                                
%  |_|  |_|\___/ \__,_|\___|_|_|_|_| |_|\__, |                                
%                                       |___/                                 

%   _    ___ ___ _ _ ____ ___                                     
%  | |  | __/ __| | |__  | __|                                    
%  | |__| _| (__|_  _|/ /|__ \                                    
%  |____|___\___| |_|/_/ |___/                                    
%                               

%   _             _       ___         _                      
%  | |   ___ _  _(_)___  / __|___  __| |_  _ _ __ _ _ _  ___ 
%  | |__/ _ \ || | (_-< | (__/ _ \/ _| ' \| '_/ _` | ' \/ -_)
%  |____\___/\_,_|_/__/  \___\___/\__|_||_|_| \__,_|_||_\___|
%                                                            

clear variables 
close all 


%% Switch to control whether model includes or exludes evaporation. 

includeEvaporation = true;

%% Paramter Setting and Setup


% Parameters for Lake A 
VwA = 40000000;      % Lake volume, m^3
FwoutA = 70000;      % Average water outflow, m^3/day (this equals the inflow)

if (includeEvaporation == true)
    FwoutA = FwoutA * (2/3); % To Account for Evaporation If Enabled.  
end


% Paramaters for Lake B 
VwB = 25000000;      % Lake volume, m^3
FwoutB = 70000;      % Average water outflow, m^3/day (this equals the inflow)

if (includeEvaporation == true)
    FwoutB = FwoutA * (3/4); % To Account for Evaporation If Enabled. 
end

% General Paramaters 
ndays = 7300;        % no. of days to simulate
Mpin = 0.15;        % Mass of pollutant entering lake per day in tonnes


%% Lake A Simulation

% Define the initial conditions, Mp0 is the initial amount of pollutant in
% the lake before the factory opens
Mp0 = 0;
% Initialise a vector Mp for storing simulated mass of pollutant in the
% lake over time
MpA = zeros(1,ndays); 
% Initialise a vector Cp for storing simulated concentration of pollutant 
% in the lake
CpA = zeros(1,ndays);  
% Calculate the first day using the initial conditions
% Mass of pollutant in lake, Mp
MpA(1) = Mp0 + Mpin - ((Mp0*FwoutA)/VwA); 
% Concentration of pollutant in lake
CpA(1) =  (MpA(1)/VwA)*10^6;
% Simulate for the remaining time period
for i=2:ndays
    MpA(i) = MpA(i-1) + Mpin - ((MpA(i-1)*FwoutA)/VwA);
    CpA(i) = (MpA(i)/VwA)*10^6;      % mass/volume *10^6 to convert to ppmv 
end


%% Lake B Simulation 

% Define the initial conditions, Mp0 is the initial amount of pollutant in
% the lake before the factory opens
Mp0 = 0;
% Initialise a vector Mp for storing simulated mass of pollutant in the
% lake over time
MpB = zeros(1,ndays); 
% Initialise a vector Cp for storing simulated concentration of pollutant 
% in the lake
CpB = zeros(1,ndays);  
% Calculate the first day using the initial conditions
% Mass of pollutant in lake, Mp
MpB(1) = Mp0 + Mpin - ((Mp0*FwoutB)/VwB); 
% Concentration of pollutant in lake
CpB(1) =  (MpB(1)/VwB)*10^6;
% Simulate for the remaining time period
for i=2:ndays
    MpB(i) = MpB(i-1) + Mpin - ((MpB(i-1)*FwoutB)/VwB);
    CpB(i) = (MpB(i)/VwB)*10^6;      % mass/volume *10^6 to convert to ppmv 
end


%% Concentration Plots

% Make a plot of concentration C, label axis, add title
% lakeAPlot = figure(1)
% plot(1:ndays, CpA)
% hold on 
% xlabel ('Number of days')
% ylabel ('Concentration (ppmv)')
% title ('Concentration of pollutant in Hyfryd Llyn Lake A ') 
% axis([0 7300 0 2.5])
% 
% lakeBPlot = figure(2)
% plot(1:ndays, CpB)
% hold on 
% xlabel ('Number of days')
% ylabel ('Concentration (ppmv)')
% title ('Concentration of pollutant in Hyfryd Llyn Lake B ') 
% axis([0 7300 0 2.5])


dualPlot = figure(3);
plot(1:ndays, CpA,'LineWidth',2)
hold on 
plot(CpB,'LineWidth',2)
xlabel ('Number of Days')
ylabel ('Concentration of Pollutant In Lake (ppmv)')
% title ('Concentration of Pollutant in Two Connected Lakes ') 
axis([0 7300 0 5])
legend(["Lake A","Lake B"])
set(gcf,'color','w');


%% Limit Finding

disp(newline) 

%%
concentrationLimit = 1.5; % Set Concentration Limit Here (ppm)

%%
% Lake A 
locationIndexA = find(CpA >= concentrationLimit);
exceedenceDayA = locationIndexA(1);
disp("Pollutant concentration limit of " + num2str(concentrationLimit) + " for Lake A is reached at day " + num2str(exceedenceDayA) + ".");

% Lake B 
locationIndexB = find(CpB >= concentrationLimit);
exceedenceDayB = locationIndexB(1);
disp("Pollutant concentration limit of " + num2str(concentrationLimit) + " for Lake B is reached at day " + num2str(exceedenceDayB) + ".");

% Highlight for Lake A 
hold on
x= [exceedenceDayA exceedenceDayA 0]; 
y= [0 concentrationLimit concentrationLimit]; %produces three xy coordinates to be plotted 
line(x,y,'Color','blue')
% text(exceedenceDayA,concentrationLimit,['\leftarrow Concentration limit is expected to be exceeded by t=' num2str(exceedenceDayA)])

% Highlight for Lake B
hold on
x= [exceedenceDayB exceedenceDayB 0]; 
y= [0 concentrationLimit concentrationLimit]; %produces three xy coordinates to be plotted 
line(x,y,'Color','red')
% text(exceedenceDayB,concentrationLimit-0.1,['\leftarrow Concentration limit is expected to be exceeded by t=' num2str(exceedenceDayB)])



%% New Concentration Input Determination

%% Set Scale Factor Here

scaleFactor = 0.8;

%%
Mp_okA = scaleFactor * (concentrationLimit / 1e6) * VwA ;
TwA = VwA/FwoutA; % the residence time of water in model
Mpin_okA = Mp_okA./TwA;    % New acceptable rate of input of pollutant to the lake

Mp_okB = scaleFactor* (concentrationLimit / 1e6) * VwB ;
TwB = VwB/FwoutB; % the residence time of water in model
Mpin_okB = Mp_okB./TwB;    % New acceptable rate of input of pollutant to the lake 

% We need to take the smallest of the acceptable inputs, to insure that the
% levels remain at the desired limit for BOTH lakes. 

if Mpin_okA < Mpin_okB
    
    Mpin_ok = Mpin_okA;
    
else
    
    Mpin_ok = Mpin_okB;
    
end


disp("To Keep Pollution " + num2str((1-scaleFactor) * 100) + "% below limit " + num2str(concentrationLimit) + ",for BOTH lakes, rate of input should be " + num2str(Mpin_ok * 10 ));


% New Simulation Lake A 

% Initialise a vector Mp_ok for storing new simulated mass of pollutant in 
% the lake.
MpA0= 0;
% Initialise a vector Mp for storing simulated mass of pollutant in the
% lake over time
MpA_ok = zeros(1,ndays); 
% Initialise a vector Cp for storing simulated concentration of pollutant 
% in the lake
CpA_ok = zeros(1,ndays);  
% Calculate the first day using the initial conditions
MpA_ok(1) = MpA0 + Mpin_ok - ((MpA0*FwoutA)/VwA); 
% Concentration of pollutant in lake
CpA_ok(1) =  (MpA_ok(1)/VwA)*10^6;
% Simulate for the time period
for i=2:ndays
    MpA_ok(i) = MpA_ok(i-1) + Mpin_ok - ((MpA_ok(i-1)*FwoutA)/VwA);
    CpA_ok(i) = (MpA_ok(i)/VwA)*10^6;      % mass/volume *10^6 to convert to ppmv 
end 

% Lake B 

% Initialise a vector Mp_ok for storing new simulated mass of pollutant in 
% the lake.
MpB0= 0;
% Initialise a vector Mp for storing simulated mass of pollutant in the
% lake over time
MpB_ok = zeros(1,ndays); 
% Initialise a vector Cp for storing simulated concentration of pollutant 
% in the lake
CpB_ok = zeros(1,ndays);  
% Calculate the first day using the initial conditions
MpB_ok(1) = MpB0 + Mpin_ok - ((MpB0*FwoutB)/VwB); 
% Concentration of pollutant in lake
CpB_ok(1) =  (MpB_ok(1)/VwB)*10^6;
% Simulate for the time period. 
for i=2:ndays
    MpB_ok(i) = MpB_ok(i-1) + Mpin_ok - ((MpB_ok(i-1)*FwoutB)/VwB);
    CpB_ok(i) = (MpB_ok(i)/VwB)*10^6;      % mass/volume *10^6 to convert to ppmv 
end 


% Plotting New Results
newConcentrationGraph = figure(4);
plot(1:ndays,CpA_ok,"blue",'LineWidth',2)
hold on 
plot(1:ndays,CpB_ok,"red",'LineWidth',2)
xlabel ('Number of Days')
ylabel ('Concentration of Pollutant in Lake (ppmv)')
% title ('Modified Concentration') 
axis([0 7300 0 1.5])
legend(["Lake A","LakeB"])
set(gcf,'color','w');


%% Extra Plotting (Currently Not Needed)


% hold on 
% plot(1:ndays,Cp_ok,"green",'LineWidth',2)
% 
% legend(["Current Rate of Pollution","New Rate of Pollution"])
% grid()
% 
% set(gcf,'color','w')
% set(gca,'color','w')
% 
% hold on 
% 
% line([0 8000],[MaxConc,MaxConc])
% text(3000,1.6,"Concentration limit")
