clear all; clc; close all;

%% -------- User Instructions --------------------------------------------
% Download this on your local drive
% Make sure you have 'importfile.m' in the same folder
% Change the directory name as you see fit
% Proceed to next section

dirN = '/Users/lykhoo/Documents/Projects/Lab/SWAPI_analysis';
cd(dirN);
species = 'H_250V';
testDate = '2024.03.27';
ISPressure = '1E-6Torr'; %For file saving
dataFolder = [testDate,'/Beam_Directionality_Test']; 
saveDir = fullfile(dirN,'Analysis');   
saveFig = 0; % Do you want to save the fig? 0 = No; 1 = Yes

%% -------- User Input: Horizontal Sweep Begining Time --------------------
% 3 keV Ar & 250 H 
% Time when moving from -ve to +ve
tBegInc     = [datetime(2024,3,27,11,45,57.7);datetime(2024,3,27,12,26,49.0);...
    datetime(2024,3,27,12,30,53.0);datetime(2024,3,27,12,42,52.3)]; 

% Time when moving from +ve to -ve
tBeg     = [datetime(2024,3,27,11,39,31.9);datetime(2024,3,27,12,23,56.5);...
    datetime(2024,3,27,12,37,28.0);datetime(2024,3,27,12,39,51.0)]; 

ctnumInc = [1,2,3,4]; % SWAPI file number that contains the horizontal sweeps. 
ctnum = [1,2,3,4]; % SWAPI file number that contains the horizontal sweeps. 


hrzCenter = 15.8244; %mm % Horizontal center for SWAPI

%% Data Analysis
SWAPICounts2022 = Arneg;
t = table2array(SWAPICounts2022(:,1));
PCEMV = table2array(SWAPICounts2022(:,2));
SCEMV = table2array(SWAPICounts2022(:,3));
ESAV  = table2array(SWAPICounts2022(:,4));
Prate = table2array(SWAPICounts2022(:,5));
Srate = table2array(SWAPICounts2022(:,6));
Crate = table2array(SWAPICounts2022(:,7));

date = datetime(t,'InputFormat','MM/dd/yyyy/HH:mm:ss.S');
horz = [-41:1:100];
hrz(1) = -41;
tStart = tBegInc(numCt);
tEnd   = tStart + seconds(length(horz));
[~,indtSt] = min(abs(date-tStart));
[~,indtEnd] = min(abs(date-tEnd));
% make the length of inner angle the same as time
trange = date(indtSt:indtEnd);
dt = seconds(trange(2:end)-trange(1:end-1));

for nt = 1:1:length(dt)
    if hrz(1) <0
        hrz(nt+1) = hrz(nt)+dt(nt);
    else
        hrz(nt+1) = hrz(nt)-dt(nt);
    end
end
pPlot = Prate(indtSt:indtEnd); 
sPlot = Srate(indtSt:indtEnd); 
cPlot = Crate(indtSt:indtEnd); 

%% --------------- Begin Analysis ---------------------------------------
for numCt = 1:length(ctnum)
    hrz = []; hrzRev = [];
    % Retrieve data
    fileN = fullfile(strcat(dirN,dataFolder), ['SWAPI-Counts_',testDate,'_',num2str(ctnumInc(numCt)),'.csv']);
    SWAPICounts2022 = importfile(fileN);

    t = table2array(SWAPICounts2022(:,1));
    PCEMV = table2array(SWAPICounts2022(:,2));
    SCEMV = table2array(SWAPICounts2022(:,3));
    ESAV  = table2array(SWAPICounts2022(:,4));
    Prate = table2array(SWAPICounts2022(:,5));
    Srate = table2array(SWAPICounts2022(:,6));
    Crate = table2array(SWAPICounts2022(:,7));

    %% Plot data from -75 to +100 mm
    date = datetime(t,'InputFormat','MM/dd/yyyy/HH:mm:ss.S');
    horz = [-41:1:100];
    hrz(1) = -41;
    tStart = tBegInc(numCt);
    tEnd   = tStart + seconds(length(horz));
    [~,indtSt] = min(abs(date-tStart));
    [~,indtEnd] = min(abs(date-tEnd));
    % make the length of inner angle the same as time
    trange = date(indtSt:indtEnd);
    dt = seconds(trange(2:end)-trange(1:end-1));

    for nt = 1:1:length(dt)
        if hrz(1) <0
            hrz(nt+1) = hrz(nt)+dt(nt);
        else
            hrz(nt+1) = hrz(nt)-dt(nt);
        end
    end
    pPlot = Prate(indtSt:indtEnd); 
    sPlot = Srate(indtSt:indtEnd); 
    cPlot = Crate(indtSt:indtEnd); 
    
    fig1 = figure(13); fig1.Position =[43 269 1443 595]; %for 5 %[193 276 1263 588]; %for 8 %[193 549 1232 315]; % for 3
    subplot(2,4,numCt)
    plot(hrz-hrzCenter,pPlot,'-','Color',[0 0.4470 0.7410],'LineWidth',2); hold on;
    plot(hrz-hrzCenter,sPlot,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2); 
    plot(hrz-hrzCenter,cPlot,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',2); 
    plot([-40 40],[1000 1000],'k:','LineWidth',1.5);
    plot([-40 40],[500 500],'k:','LineWidth',1.5);
    plot([-40 40],[100 100],'k:','LineWidth',1.5);
    plot([0 0],[100 100000],'k--','LineWidth',1.5);
    xlabel('Horizontal (mm)'); ylabel('Counts (Hz)'); 
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',12,'TickLength',[0.02 0.1],'YScale','log');
    leg = legend('PCEM','SCEM','Coin','NumColumns',3); leg.Box = 'off'; grid on; grid minor;
    leg.Location = 'northwest'; grid on;
    switch numCt
        case 1 
        title('Helmholtz Off: 3keV Argon');
        case 2
        title('Helmholtz On: 0.60A');
        case 3
        title('Helmholtz On: 0.61A');
        case 4
        title('Helmholtz On: 0.62A');
        case 5
        title('Helmholtz On: 0.63A');
        case 6
        title('Helmholtz On: 0.64A');
        case 7
        title('Helmholtz On: 0.65A');
        case 8
        title('Helmholtz On: 0.66A');
    end
    ylim([100 2000]); xlim([-40 40]);


%% Plot data from +100 to -75 mm
    % Retrieve data
    fileN = fullfile(strcat(dirN,dataFolder), ['SWAPI-Counts_',testDate,'_',num2str(ctnum(numCt)),'.csv']);
    SWAPICounts2022 = importfile_v2(fileN);

    t = table2array(SWAPICounts2022(:,1));
    PCEMV = table2array(SWAPICounts2022(:,2));
    SCEMV = table2array(SWAPICounts2022(:,3));
    ESAV  = table2array(SWAPICounts2022(:,4));
    Prate = table2array(SWAPICounts2022(:,5));
    Srate = table2array(SWAPICounts2022(:,6));
    Crate = table2array(SWAPICounts2022(:,7));

    date = datetime(t,'InputFormat','MM/dd/yyyy/HH:mm:ss.S');
    horz = [-75:1:100];
    hrzRev(1) = 100;
    tStart = tBeg(numCt);
    tEnd   = tStart + seconds(length(horz));
    [~,indtSt] = min(abs(date-tStart));
    [~,indtEnd] = min(abs(date-tEnd));
    % make the length of inner angle the same as time
    trange = date(indtSt:indtEnd);
    dt = seconds(trange(2:end)-trange(1:end-1));

    for nt = 1:1:length(dt)
        if hrzRev(1) <0
            hrzRev(nt+1) = hrzRev(nt)+dt(nt);
        else
            hrzRev(nt+1) = hrzRev(nt)-dt(nt);
        end
    end
    pPlot = Prate(indtSt:indtEnd); 
    sPlot = Srate(indtSt:indtEnd); 
    cPlot = Crate(indtSt:indtEnd); 
    
    fig2 = figure(15); fig2.Position =[43 269 1443 595]; %for 5 %[193 276 1263 588]; %for 8 %[193 549 1232 315]; % for 3
    subplot(2,4,numCt)
    plot(hrzRev-hrzCenter,pPlot,'-','Color',[0 0.4470 0.7410],'LineWidth',2); hold on;
    plot(hrzRev-hrzCenter,sPlot,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2); 
    plot(hrzRev-hrzCenter,cPlot,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',2); 
    plot([-40 40],[1000 1000],'k:','LineWidth',1.5);
    plot([-40 40],[500 500],'k:','LineWidth',1.5);
    plot([-40 40],[100 100],'k:','LineWidth',1.5);
    plot([0 0],[100 100000],'k--','LineWidth',1.5);
    xlabel('Horizontal (mm)'); ylabel('Counts (Hz)'); 
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',12,'TickLength',[0.02 0.1],'YScale','log');
    leg = legend('PCEM','SCEM','Coin','NumColumns',3); leg.Box = 'off'; grid on; grid minor;
    leg.Location = 'northwest'; grid on;
    switch numCt
        case 1 
        title('Helmholtz Off: 3keV Argon');
        case 2
        title('Helmholtz On: 0.60A');
        case 3
        title('Helmholtz On: 0.61A');
        case 4
        title('Helmholtz On: 0.62A');
        case 5
        title('Helmholtz On: 0.63A');
        case 6
        title('Helmholtz On: 0.64A');
        case 7
        title('Helmholtz On: 0.65A');
        case 8
        title('Helmholtz On: 0.66A');
    end
    ylim([100 2000]); xlim([-40 40]);
end
%% save figures
if saveFig == 1
    fileName1 = ['HorizontalSweep__fnum',num2str(ctnum(1)),'_',num2str(ctnum(end)),'_all_Horizontal_-75to100mm_1kV_',species,'_',ISPressure,testDate];
    saveFullDir = fullfile(saveDir,filename1);
    saveas(fig1,saveFullDir,'png');

    filename2 = ['HorizontalSweep__fnum',num2str(ctnum(1)),'_',num2str(ctnum(end)),'_all_Horizontal_100to-75mm_1kV_',species,'_',ISPressure,testDate];
    saveFullDir = fullfile(saveDir,filename2);
    saveas(fig2,saveFullDir,'png');
end