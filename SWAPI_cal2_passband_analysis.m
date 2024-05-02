clear all; clc; format long;
format long

%% Initial setup

% ===== Controlled variables ===========================
beam_energy = 1.0; % unit: kV
k_factor    = 1.88; % ref : SWAP data
dotsize     = 60;   % For passband purposes. nominal = 60 & extended = 20
save_plot   = 0;    % Yes:1 or No:0, default: 0
save_data   = 0;    % Yes:1 or No:0, default: 0
disp_2x3set = 2;    % PRM/SEC/COI and flip Az : 1
                    % PRM/SEC/COI and ADE, S/P: 2
                    % If SWAPIcdr = 1, we don't care this "disp_set" value
advanced    = 1;    % Yes:1 or No:0, default: 0
disp_1x3set = 1;    % 1: COI/ADE/SPR ; 2: PRM/SEC/COI

% ===== Set directories ================================
dir = '/Users/lykhoo/Library/CloudStorage/OneDrive-PrincetonUniversity/Documents/Lab/SWAPI Cal/FM-Cal/';
save_folder = [dir '/Analysis/'];
% fld = ['FM_Pre_CAL_databaseCal/Raw_Data/'];
fld = 'Data_analysis/';

% ===== Set X-ticks and Y-ticks values =================
ESA_V_range = (-420:-10:-640) * (beam_energy);
% ESA_V_range = [-750:-30:-2400]; % for extended ESA sweep
OUTER_range = [0 -15:1:15 0];

% ===== Import database ================================
tdate = '2022.05.10';
dd = str2num(tdate(end-1:end)); mm = str2num(tdate(end-3:end-2));
yyyy = str2num(tdate(1:4));
passbandFolder    = 'passband-1_He_sunglass/';
name_range  = 1:1:33; % Passband Name Range

clc; disp('Initial Setup Completed')


%% Import data and binding
% ===== Set readtable format ===========================
opts = detectImportOptions([dir fld passbandFolder 'SWEEP_ESA_',tdate,'_',num2str(name_range(1)),'.csv']);
opts.DataLines = 3; opts.VariableNamesLine = 2;

for i = 1:1:length(name_range)
    name_num = name_range(i);
    data = readtable([dir fld passbandFolder 'SWEEP_ESA_',tdate,'_',num2str(name_range(i)),'.csv'] ...
        ,opts,'ReadVariableNames',false);
    ESA_V_str(:,i) = string(table2array(data(2:end-2,3)));
    for j = 1:size(ESA_V_str,1)
        tmp = char(ESA_V_str(j,i));
        ESA_V_bin(j,i) = str2num(tmp(6:9));
        clear tmp;
    end

    I_PRM_bin(:,i) = table2array(data(2:end-2,8));
    I_SEC_bin(:,i) = table2array(data(2:end-2,9));
    I_COI_bin(:,i) = table2array(data(2:end-2,10));
    I_S2P_bin(:,i) = table2array(data(2:end-2,9))./table2array(data(2:end-2,8));
    I_P2S_bin(:,i) = table2array(data(2:end-2,8))./table2array(data(2:end-2,9));
    
    clear name_num data
end

% ===== Data binning & processing ======================
% I_PRM_bin(isnan(I_PRM_bin)) = 0;
% I_SEC_bin(isnan(I_SEC_bin)) = 0;
% I_COI_bin(isnan(I_COI_bin)) = 0;
% I_S2P_bin(isnan(I_S2P_bin)) = 0;
% I_S2P(isinf(I_S2P)) = NaN;

I_PRM_norm = I_PRM_bin./max(max(I_PRM_bin));
I_SEC_norm = I_SEC_bin./max(max(I_SEC_bin));
I_COI_norm = I_COI_bin./max(max(I_COI_bin));

I_PRM_norm_log = log10(I_PRM_bin)./max(max(log10(I_PRM_bin)));
I_SEC_norm_log = log10(I_SEC_bin)./max(max(log10(I_SEC_bin)));
I_COI_norm_log = log10(I_COI_bin)./max(max(log10(I_COI_bin)));

% I_PRM_norm_log(isinf(I_PRM_norm_log)) = NaN;
% I_SEC_norm_log(isinf(I_SEC_norm_log)) = NaN;
% I_COI_norm_log(isinf(I_COI_norm_log)) = NaN;

% ===== For scatter plot: data binding =================
% ESA voltage
ESA_V_range_matrix = ones(1,size(I_PRM_norm_log,2)) .* ESA_V_range';
ESA_V_range_collapse = reshape(ESA_V_range_matrix,1,[]);

% outer rotational angle
OUTER_range_matrix = ones(1,size(I_PRM_norm_log,1)) .* OUTER_range';
OUTER_range_collapse = reshape(OUTER_range_matrix',1,[]);

% PRM/SEC/COI rate
I_PRM_norm_collapse = reshape(I_PRM_norm,1,[]);
I_PRM_norm_log_collapse = reshape(I_PRM_norm_log,1,[]);
I_PRM_collapse = reshape(I_PRM_bin,1,[]);
I_SEC_norm_collapse = reshape(I_SEC_norm,1,[]);
I_SEC_norm_log_collapse = reshape(I_SEC_norm_log,1,[]);
I_SEC_collapse = reshape(I_SEC_bin,1,[]);
I_COI_norm_collapse = reshape(I_COI_norm,1,[]);
I_COI_norm_log_collapse = reshape(I_COI_norm_log,1,[]);
I_COI_collapse = reshape(I_COI_bin,1,[]);

% advanced calculations
abs_detect_eff = I_COI_collapse.^2./(I_PRM_collapse.*I_SEC_collapse);
PRM_SEC_ratio  = I_PRM_collapse./I_SEC_collapse;
SEC_PRM_ratio  = I_SEC_collapse./I_PRM_collapse;

% output saved .mat file if needed
if save_data == 1
    save(['SWAPI_Cal2_inner_' num2str(inner_angle,'%02d') '.mat'],...
            'I_PRM_bin','I_SEC_bin','I_COI_bin','I_S2P_bin','I_P2S_bin',...
            'OUTER_range','ESA_V_range',...
            'I_PRM_norm','I_SEC_norm','I_COI_norm',...
            'I_PRM_norm_log','I_SEC_norm_log','I_COI_norm_log',...
            'ESA_V_range_collapse','OUTER_range_collapse',...
            'I_PRM_norm_log_collapse','I_PRM_collapse',...
            'I_SEC_norm_log_collapse','I_SEC_collapse',...
            'I_COI_norm_log_collapse','I_COI_collapse',...
            'abs_detect_eff','PRM_SEC_ratio','SEC_PRM_ratio')
else
    
end
 

%% Norm plot (upper:raw; lower:process)
close all; clc;

% ===== "normal" 1x3 angle/ESA raw on top, 1x3 angle/beamE raw at bottom
if advanced ~= 1 

    % ===== Set figure parameters ========
    f2 = figure(2);
    set(gcf,'position',[20,100,1000,600])
    f2.Color = 'w';

    % ===== PCEM raw rate normalized =====
    sf1 = subplot(2,3,1);
    sf1.Position = [0.06 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,I_PRM_norm_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on; 
    ylabel(h, 'Normalized Counts','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Primary CEM counts (raw)')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== SCEM raw rate normalized =====
    sf2 = subplot(2,3,2);
    sf2.Position = [0.39 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,I_SEC_norm_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Counts','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary CEM counts (raw)')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== COIN raw rate normalized =====
    sf3 = subplot(2,3,3);
    sf3.Position = [0.72 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,I_COI_norm_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Counts','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence counts (raw)')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);


    switch disp_2x3set
        case 1 % raw bottom with beam energy
    
    % ===== PCEM raw rate normalized conversion to energy =====
    sf4 = subplot(2,3,4);
    sf4.Position = [0.06 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse*-k_factor,-OUTER_range_collapse,dotsize,I_PRM_norm_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; grid on; box on;
    ylabel(h, 'Normalized Counts','FontWeight','bold','FontSize',11)
    xlabel('Ion Energy [eV]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([770*beam_energy 1220*beam_energy]); 
    xticks([beam_energy*1e3*.8:beam_energy*40:beam_energy*1e3*1.2]);
    % ===== for extended ESA sweep =====
    % xlim([440*beam_energy 1540*beam_energy]); 
    % xticks([beam_energy*1e3*.5:beam_energy*80:beam_energy*1e3*1.5]);
    % ==================================
    % tickarray = flip(round([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]*-k_factor));
    % xticks(tickarray);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Primary CEM counts')
    % set(gca, 'XDir','reverse');
    set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    
    % ===== SCEM raw counts normalized conversion to energy =====
    sf5 = subplot(2,3,5);
    sf5.Position = [0.39 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse*-k_factor,-OUTER_range_collapse,dotsize,I_SEC_norm_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; grid on; box on;
    ylabel(h, 'Normalized Counts','FontWeight','bold','FontSize',11)
    xlabel('Ion Energy [eV]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([770*beam_energy 1220*beam_energy]); 
    xticks([beam_energy*1e3*.8:beam_energy*40:beam_energy*1e3*1.2]);
    % ===== for extended ESA sweep =====
    % xlim([440*beam_energy 1540*beam_energy]);
    % xticks([beam_energy*1e3*.5:beam_energy*80:beam_energy*1e3*1.5]);
    % ==================================
    % tickarray = flip(round([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]*-k_factor));
    % xticks(tickarray);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary CEM counts')
    % set(gca, 'XDir','reverse'); 
    set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    
    % ===== COIN raw counts normalized conversion to energy =====
    sf6 = subplot(2,3,6);
    sf6.Position = [0.72 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse*-k_factor,-OUTER_range_collapse,dotsize,I_COI_norm_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; grid on; box on;
    ylabel(h, 'Normalized Counts','FontWeight','bold','FontSize',11)
    xlabel('Ion Energy [eV]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([770*beam_energy 1220*beam_energy]); 
    xticks([beam_energy*1e3*.8:beam_energy*40:beam_energy*1e3*1.2]);
    % ===== for extended ESA sweep =====
    % xlim([440*beam_energy 1540*beam_energy]);
    % xticks([beam_energy*1e3*.5:beam_energy*80:beam_energy*1e3*1.5]);
    % ==================================
    % tickarray = flip(round([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]*-k_factor));
    % xticks(tickarray);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence counts')
    % set(gca, 'XDir','reverse'); 
    set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % =================================================================
    % ===== that have Abs Detection Efficiency and S/P ratio ==========
    % =================================================================
        case 2 
    % ===== Abs Detection Efficiency =====
    sf4 = subplot(2,3,4);
    sf4.Position = [0.06 0.12 0.25 0.34];
    idx_zero = find(abs_detect_eff==0);
    idx_one  = find(abs_detect_eff==1);
    abs_detect_eff_filtered = abs_detect_eff; 
    abs_detect_eff_filtered(idx_zero) = nan; abs_detect_eff_filtered(idx_one) = nan;
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,abs_detect_eff_filtered,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_zero),OUTER_range_collapse(idx_zero),dotsize*.8,...
            abs_detect_eff(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor','g')
    scatter(ESA_V_range_collapse(idx_one),OUTER_range_collapse(idx_one),dotsize*.8,...
            abs_detect_eff(idx_one),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor','m')

    set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.02:.2); colormap turbo; box on; 
    ylabel(h, 'C^2/(P*S)','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([.01 .2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Absolute Detection Efficiency')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf);
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  

    annotation('textbox',[0.29, 0.50, 0.08, 0],'string','1'     ,'FontSize',11,'LineStyle','none');
    annotation('textbox',[0.29, 0.11, 0.08, 0],'string','No COI','FontSize',11,'LineStyle','none');
    annotation('rectangle',[0.28, 0.48 .006 .01],'FaceColor','m')  
    annotation('rectangle',[0.28, 0.09 .006 .01],'FaceColor','g') 

    
    % ===== S/P ratio =====
    sf5 = subplot(2,3,5);
    sf5.Position = [0.39 0.12 0.25 0.34];
    idx_zero = find(SEC_PRM_ratio==0);
    idx_inf  = find(SEC_PRM_ratio==inf);
    SEC_PRM_ratio_filtered = SEC_PRM_ratio; 
    SEC_PRM_ratio_filtered(idx_zero) = nan; SEC_PRM_ratio_filtered(idx_inf) = nan;
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,SEC_PRM_ratio_filtered,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_zero),OUTER_range_collapse(idx_zero),dotsize*.8,...
            SEC_PRM_ratio(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor','g')
    scatter(ESA_V_range_collapse(idx_inf),OUTER_range_collapse(idx_inf),dotsize*.8,...
            SEC_PRM_ratio(idx_inf),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.4940 0.1840 0.5560])

    % set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.2:2); colormap(sf5,turbo); box on; 
    ylabel(h, 'S/P Ratio','FontWeight','bold','FontSize',11)
    % set(gca,'ColorScale','log')
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary/Primary Ratio')
    set(gca, 'XDir','reverse'); %set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf); grid on;
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  

    annotation('textbox',[0.62, 0.50, 0.08, 0],'string','SEC only','FontSize',11,'LineStyle','none'); 
    annotation('textbox',[0.62, 0.11, 0.08, 0],'string','PRM only','FontSize',11,'LineStyle','none'); 
    annotation('rectangle',[0.61, 0.48 .006 .01],'FaceColor',[0.4940 0.1840 0.5560])  
    annotation('rectangle',[0.61, 0.09 .006 .01],'FaceColor','g')  
    
    end


elseif advanced == 1 % For 1x3 display set, SWAPI CDR slides

    switch disp_1x3set
        case 1 % COI/ADE/SPR
    
    f3 = figure(3);
    set(gcf,'position',[20,100,1000,300])
    f3.Color = 'w';

    % ===== COIN raw rate normalized =====
    sf3 = subplot(1,3,1);
    sf3.Position = [0.05 0.22 0.25 0.68]; %[0.72 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_COI_norm_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== absolute detection efficiency =====
    sf4 = subplot(1,3,2);
    sf4.Position = [0.38 0.22 0.25 0.68];
    idx_PSonly  = find(I_COI_collapse==0 & (I_SEC_collapse ~=0 | I_PRM_collapse ~=0));
    idx_zero    = find(abs_detect_eff==0);
    idx_one     = find(abs_detect_eff==1);
    abs_detect_eff_filtered = abs_detect_eff; 
    abs_detect_eff_filtered(idx_zero) = nan; abs_detect_eff_filtered(idx_one) = nan;
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,abs_detect_eff_filtered,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_PSonly),-OUTER_range_collapse(idx_PSonly),dotsize*.8,...
            abs_detect_eff(idx_PSonly),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','x','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',.8)
    scatter(ESA_V_range_collapse(idx_zero),-OUTER_range_collapse(idx_zero),dotsize*.8,...
            abs_detect_eff(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','x','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',.8)
    scatter(ESA_V_range_collapse(idx_one),-OUTER_range_collapse(idx_one),dotsize*.8,...
            abs_detect_eff(idx_one),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','+','MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',.8)

    set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.02:.2); colormap turbo; box on; 
    ylabel(h, 'C^2/(P*S)','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([.01 .2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Absolute Detection Efficiency')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf);
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  
    grid on; sf4.GridAlpha = .05;

    annotation('textbox', [0.608, 0.98, 0, 0],'string','1','FontSize',11);
    % annotation('rectangle',[0.60, 0.94 .005 .018],'FaceColor','m')  
    annotation('textbox', [0.608, 0.21, 0.08, 0],'string','No COI','FontSize',11,'LineStyle','none');
    % annotation('rectangle',[0.60, 0.17 .005 .018],'FaceColor','b') 

    ax = gca; % get the current axis
    ax.Clipping = 'off';
    scatter(min(ESA_V_range)-40*beam_energy,-18,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','x',...
            'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1.5)
    scatter(min(ESA_V_range)-40*beam_energy,+18.2,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','+',...
            'MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',1.5)


    % ===== SEC/PRM ratio =====
    sf5 = subplot(1,3,3);
    sf5.Position = [0.71 0.22 0.25 0.68]; %[0.39 0.12 0.25 0.34];
    idx_zero = find(SEC_PRM_ratio==0);
    idx_inf  = find(SEC_PRM_ratio==inf);
    SEC_PRM_ratio_filtered = SEC_PRM_ratio; 
    SEC_PRM_ratio_filtered(idx_zero) = nan; SEC_PRM_ratio_filtered(idx_inf) = nan;
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,SEC_PRM_ratio_filtered,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_zero),-OUTER_range_collapse(idx_zero),dotsize*.5,...
            SEC_PRM_ratio(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,...
            'marker','v','MarkerEdgeColor','b','MarkerFaceColor','none','LineWidth',.7)
    scatter(ESA_V_range_collapse(idx_inf),-OUTER_range_collapse(idx_inf),dotsize*.5,...
            SEC_PRM_ratio(idx_inf),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,...
            'marker','^','MarkerEdgeColor','m','MarkerFaceColor','none','LineWidth',.7)

    % set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.2:2); colormap(sf5,turbo); box on; 
    ylabel(h, 'S/P Ratio','FontSize',12)
    % set(gca,'ColorScale','log')
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([.1 2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary/Primary Ratio')
    set(gca, 'XDir','reverse'); %set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf); grid on;
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  
    grid on; sf5.GridAlpha = .05;

    annotation('textbox',[0.9385, 0.98, 0.08, 0],'string','No PRM','FontSize',11,'LineStyle','none'); 
    % annotation('rectangle',[0.93, 0.94 .005 .018],'FaceColor','m')
    annotation('textbox',[0.9385, 0.21, 0.08, 0],'string','No SEC','FontSize',11,'LineStyle','none'); 
    % annotation('rectangle',[0.93, 0.17 .005 .018],'FaceColor','b')

    ax = gca; % get the current axis
    ax.Clipping = 'off';
    scatter(min(ESA_V_range)-40*beam_energy,-17.8,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,'marker','v',...
            'MarkerEdgeColor','b','MarkerFaceColor','none','LineWidth',.7)
    scatter(min(ESA_V_range)-40*beam_energy,+18.2,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,'marker','^',...
            'MarkerEdgeColor','m','MarkerFaceColor','none','LineWidth',.7)

        
        case 2 % PRM/SEC/COI for CDR (Az angle vs ESA)

    f4 = figure(4);
    set(gcf,'position',[20,100,1000,300])
    f4.Color = 'w';

    % ===== PRM rate normalized =====
    sf1 = subplot(1,3,1);
    sf1.Position = [0.05 0.22 0.25 0.68]; %[0.72 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_PRM_norm_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Primary Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== SEC rate normalized =====
    sf2 = subplot(1,3,2);
    sf2.Position = [0.38 0.22 0.25 0.68];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_SEC_norm_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);    

    % ===== COI rate normalized =====
    sf3 = subplot(1,3,3);
    sf3.Position = [0.71 0.22 0.25 0.68]; %[0.39 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_COI_norm_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);      
  
    end
    
end

% ===== Parameters for PLOT SAVING ===== 
set(0, 'DefaultFigurePaperPositionMode', 'auto')
set(gcf, 'InvertHardcopy', 'off')
if save_plot == 1
    print('-dpng'  ,'-r300', ['/Users/lykhoo/Library/CloudStorage/OneDrive-PrincetonUniversity/Cal4.0_database/Data Analysis/' ...
        datefld(11:14) '-' datefld(16:17) '-' num2str(dd,'%02d') ...
          '_halfPassband_InnRot_' num2str(inner_angle) '_500eV_Ar_5E-8torr_old_IS_parameters.png'])
else
end
disp('Mission Completed')

%% Log-norm plot (upper:raw; lower:process)
close all; clc;

% ===== "normal" 1x3 angle/ESA raw on top, 1x3 angle/beamE raw at bottom
if advanced ~= 1 

    % ===== Set figure parameters ========
    f2 = figure(2);
    set(gcf,'position',[20,100,1000,600])
    f2.Color = 'w';

    % ===== PCEM raw rate normalized =====
    sf1 = subplot(2,3,1);
    sf1.Position = [0.06 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,I_PRM_norm_log_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on; 
    ylabel(h, 'Normalized Log Counts','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Primary CEM counts (raw)')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== SCEM raw rate normalized =====
    sf2 = subplot(2,3,2);
    sf2.Position = [0.39 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,I_SEC_norm_log_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Log Counts','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary CEM counts (raw)')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== COIN raw rate normalized =====
    sf3 = subplot(2,3,3);
    sf3.Position = [0.72 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,I_COI_norm_log_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Log Counts','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence counts (raw)')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);


    switch disp_2x3set
        case 1 % raw bottom with beam energy
    
    % ===== PCEM raw rate normalized conversion to energy =====
    sf4 = subplot(2,3,4);
    sf4.Position = [0.06 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse*-k_factor,-OUTER_range_collapse,dotsize,I_PRM_norm_log_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; grid on; box on;
    ylabel(h, 'Normalized Log Counts','FontWeight','bold','FontSize',11)
    xlabel('Ion Energy [eV]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([770*beam_energy 1220*beam_energy]); 
    xticks([beam_energy*1e3*.8:beam_energy*40:beam_energy*1e3*1.2]);
    % ===== for extended ESA sweep =====
    % xlim([440*beam_energy 1540*beam_energy]); 
    % xticks([beam_energy*1e3*.5:beam_energy*80:beam_energy*1e3*1.5]);
    % ==================================
    % tickarray = flip(round([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]*-k_factor));
    % xticks(tickarray);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Primary CEM counts')
    % set(gca, 'XDir','reverse');
    set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    
    % ===== SCEM raw counts normalized conversion to energy =====
    sf5 = subplot(2,3,5);
    sf5.Position = [0.39 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse*-k_factor,-OUTER_range_collapse,dotsize,I_SEC_norm_log_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; grid on; box on;
    ylabel(h, 'Normalized Log Counts','FontWeight','bold','FontSize',11)
    xlabel('Ion Energy [eV]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([770*beam_energy 1220*beam_energy]); 
    xticks([beam_energy*1e3*.8:beam_energy*40:beam_energy*1e3*1.2]);
    % ===== for extended ESA sweep =====
    % xlim([440*beam_energy 1540*beam_energy]);
    % xticks([beam_energy*1e3*.5:beam_energy*80:beam_energy*1e3*1.5]);
    % ==================================
    % tickarray = flip(round([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]*-k_factor));
    % xticks(tickarray);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary CEM counts')
    % set(gca, 'XDir','reverse'); 
    set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    
    % ===== COIN raw counts normalized conversion to energy =====
    sf6 = subplot(2,3,6);
    sf6.Position = [0.72 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse*-k_factor,-OUTER_range_collapse,dotsize,I_COI_norm_log_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; grid on; box on;
    ylabel(h, 'Normalized Log Counts','FontWeight','bold','FontSize',11)
    xlabel('Ion Energy [eV]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([770*beam_energy 1220*beam_energy]); 
    xticks([beam_energy*1e3*.8:beam_energy*40:beam_energy*1e3*1.2]);
    % ===== for extended ESA sweep =====
    % xlim([440*beam_energy 1540*beam_energy]);
    % xticks([beam_energy*1e3*.5:beam_energy*80:beam_energy*1e3*1.5]);
    % ==================================
    % tickarray = flip(round([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]*-k_factor));
    % xticks(tickarray);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence counts')
    % set(gca, 'XDir','reverse'); 
    set(gca, 'YDir','normal');
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % =================================================================
    % ===== that have Abs Detection Efficiency and S/P ratio ==========
    % =================================================================
        case 2 
    % ===== Abs Detection Efficiency =====
    sf4 = subplot(2,3,4);
    sf4.Position = [0.06 0.12 0.25 0.34];
    idx_zero = find(abs_detect_eff==0);
    idx_one  = find(abs_detect_eff==1);
    abs_detect_eff_filtered = abs_detect_eff; 
    abs_detect_eff_filtered(idx_zero) = nan; abs_detect_eff_filtered(idx_one) = nan;
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,abs_detect_eff_filtered,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_zero),OUTER_range_collapse(idx_zero),dotsize*.8,...
            abs_detect_eff(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor','g')
    scatter(ESA_V_range_collapse(idx_one),OUTER_range_collapse(idx_one),dotsize*.8,...
            abs_detect_eff(idx_one),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor','m')

    set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.02:.2); colormap turbo; box on; 
    ylabel(h, 'C^2/(P*S)','FontWeight','bold','FontSize',11)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([.01 .2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Absolute Detection Efficiency')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf);
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  

    annotation('textbox',[0.29, 0.50, 0.08, 0],'string','1'     ,'FontSize',11,'LineStyle','none');
    annotation('textbox',[0.29, 0.11, 0.08, 0],'string','No COI','FontSize',11,'LineStyle','none');
    annotation('rectangle',[0.28, 0.48 .006 .01],'FaceColor','m')  
    annotation('rectangle',[0.28, 0.09 .006 .01],'FaceColor','g') 

    
    % ===== S/P ratio =====
    sf5 = subplot(2,3,5);
    sf5.Position = [0.39 0.12 0.25 0.34];
    idx_zero = find(SEC_PRM_ratio==0);
    idx_inf  = find(SEC_PRM_ratio==inf);
    SEC_PRM_ratio_filtered = SEC_PRM_ratio; 
    SEC_PRM_ratio_filtered(idx_zero) = nan; SEC_PRM_ratio_filtered(idx_inf) = nan;
    scatter(ESA_V_range_collapse,OUTER_range_collapse,dotsize,SEC_PRM_ratio_filtered,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_zero),OUTER_range_collapse(idx_zero),dotsize*.8,...
            SEC_PRM_ratio(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor','g')
    scatter(ESA_V_range_collapse(idx_inf),OUTER_range_collapse(idx_inf),dotsize*.8,...
            SEC_PRM_ratio(idx_inf),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','s','MarkerEdgeColor','k','MarkerFaceColor',[0.4940 0.1840 0.5560])

    % set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.2:2); colormap(sf5,turbo); box on; 
    ylabel(h, 'S/P Ratio','FontWeight','bold','FontSize',11)
    % set(gca,'ColorScale','log')
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Outer Rotational Stage [deg]','FontWeight','bold')
    caxis([0 2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary/Primary Ratio')
    set(gca, 'XDir','reverse'); %set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf); grid on;
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  

    annotation('textbox',[0.62, 0.50, 0.08, 0],'string','SEC only','FontSize',11,'LineStyle','none'); 
    annotation('textbox',[0.62, 0.11, 0.08, 0],'string','PRM only','FontSize',11,'LineStyle','none'); 
    annotation('rectangle',[0.61, 0.48 .006 .01],'FaceColor',[0.4940 0.1840 0.5560])  
    annotation('rectangle',[0.61, 0.09 .006 .01],'FaceColor','g')  
    
    end


elseif advanced == 1 % For 1x3 display set, SWAPI CDR slides

    switch disp_1x3set
        case 1 % COI/ADE/SPR
    
    f3 = figure(3);
    set(gcf,'position',[20,100,1000,300])
    f3.Color = 'w';

    % ===== COIN raw rate normalized =====
    sf3 = subplot(1,3,1);
    sf3.Position = [0.05 0.22 0.25 0.68]; %[0.72 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_COI_norm_log_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Log Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== absolute detection efficiency =====
    sf4 = subplot(1,3,2);
    sf4.Position = [0.38 0.22 0.25 0.68];
    idx_PSonly  = find(I_COI_collapse==0 & (I_SEC_collapse ~=0 | I_PRM_collapse ~=0));
    idx_zero    = find(abs_detect_eff==0);
    idx_one     = find(abs_detect_eff==1);
    abs_detect_eff_filtered = abs_detect_eff; 
    abs_detect_eff_filtered(idx_zero) = nan; abs_detect_eff_filtered(idx_one) = nan;
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,abs_detect_eff_filtered,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_PSonly),-OUTER_range_collapse(idx_PSonly),dotsize*.8,...
            abs_detect_eff(idx_PSonly),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','x','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',.8)
    scatter(ESA_V_range_collapse(idx_zero),-OUTER_range_collapse(idx_zero),dotsize*.8,...
            abs_detect_eff(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','x','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',.8)
    scatter(ESA_V_range_collapse(idx_one),-OUTER_range_collapse(idx_one),dotsize*.8,...
            abs_detect_eff(idx_one),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'marker','+','MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',.8)

    set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.02:.2); colormap turbo; box on; 
    ylabel(h, 'C^2/(P*S)','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([.01 .2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Absolute Detection Efficiency')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf);
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  
    grid on; sf4.GridAlpha = .05;

    annotation('textbox', [0.608, 0.98, 0, 0],'string','1','FontSize',11);
    % annotation('rectangle',[0.60, 0.94 .005 .018],'FaceColor','m')  
    annotation('textbox', [0.608, 0.21, 0.08, 0],'string','No COI','FontSize',11,'LineStyle','none');
    % annotation('rectangle',[0.60, 0.17 .005 .018],'FaceColor','b') 

    ax = gca; % get the current axis
    ax.Clipping = 'off';
    scatter(min(ESA_V_range)-40*beam_energy,-18,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','x',...
            'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1.5)
    scatter(min(ESA_V_range)-40*beam_energy,+18.2,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','+',...
            'MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',1.5)


    % ===== SEC/PRM ratio =====
    sf5 = subplot(1,3,3);
    sf5.Position = [0.71 0.22 0.25 0.68]; %[0.39 0.12 0.25 0.34];
    idx_zero = find(SEC_PRM_ratio==0);
    idx_inf  = find(SEC_PRM_ratio==inf);
    SEC_PRM_ratio_filtered = SEC_PRM_ratio; 
    SEC_PRM_ratio_filtered(idx_zero) = nan; SEC_PRM_ratio_filtered(idx_inf) = nan;
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,SEC_PRM_ratio_filtered,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o'); hold on;
    scatter(ESA_V_range_collapse(idx_zero),-OUTER_range_collapse(idx_zero),dotsize*.5,...
            SEC_PRM_ratio(idx_zero),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,...
            'marker','v','MarkerEdgeColor','b','MarkerFaceColor','none','LineWidth',.7)
    scatter(ESA_V_range_collapse(idx_inf),-OUTER_range_collapse(idx_inf),dotsize*.5,...
            SEC_PRM_ratio(idx_inf),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,...
            'marker','^','MarkerEdgeColor','m','MarkerFaceColor','none','LineWidth',.7)

    % set(gca,'color',1*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.2:2); colormap(sf5,turbo); box on; 
    ylabel(h, 'S/P Ratio','FontSize',12)
    % set(gca,'ColorScale','log')
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([.1 2])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary/Primary Ratio')
    set(gca, 'XDir','reverse'); %set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12); 
    % a = gca; b = copyobj(a, gcf); grid on;
    % set(b,'Xcolor',0*[1 1 1],'YColor',0*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
    %       'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);  
    grid on; sf5.GridAlpha = .05;

    annotation('textbox',[0.9385, 0.98, 0.08, 0],'string','No PRM','FontSize',11,'LineStyle','none'); 
    % annotation('rectangle',[0.93, 0.94 .005 .018],'FaceColor','m')
    annotation('textbox',[0.9385, 0.21, 0.08, 0],'string','No SEC','FontSize',11,'LineStyle','none'); 
    % annotation('rectangle',[0.93, 0.17 .005 .018],'FaceColor','b')

    ax = gca; % get the current axis
    ax.Clipping = 'off';
    scatter(min(ESA_V_range)-40*beam_energy,-17.8,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,'marker','v',...
            'MarkerEdgeColor','b','MarkerFaceColor','none','LineWidth',.7)
    scatter(min(ESA_V_range)-40*beam_energy,+18.2,dotsize*.5,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',.6,'marker','^',...
            'MarkerEdgeColor','m','MarkerFaceColor','none','LineWidth',.7)

        
        case 2 % PRM/SEC/COI for CDR (Az angle vs ESA)

    f4 = figure(4);
    set(gcf,'position',[20,100,1000,300])
    f4.Color = 'w';

    % ===== PRM rate normalized =====
    sf1 = subplot(1,3,1);
    sf1.Position = [0.05 0.22 0.25 0.68]; %[0.72 0.6 0.25 0.34];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_PRM_norm_log_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Log Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Primary Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

    % ===== SEC rate normalized =====
    sf2 = subplot(1,3,2);
    sf2.Position = [0.38 0.22 0.25 0.68];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_SEC_norm_log_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Log Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Secondary Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);    

    % ===== COI rate normalized =====
    sf3 = subplot(1,3,3);
    sf3.Position = [0.71 0.22 0.25 0.68]; %[0.39 0.12 0.25 0.34];
    scatter(ESA_V_range_collapse,-OUTER_range_collapse,dotsize,I_COI_norm_log_collapse,'filled',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','o')
    set(gca,'color',0*[1 1 1]); hold on; 
    h = colorbar('XTick', 0:.1:1); colormap turbo; box on;
    ylabel(h, 'Normalized Log Rate','FontSize',12)
    xlabel('ESA Voltage [V]','FontWeight','bold')
    ylabel('Azimuthal Angle [deg]','FontWeight','bold')
    caxis([0 1])
    xlim([min(ESA_V_range)-20*beam_energy max(ESA_V_range)+20*beam_energy]); 
    xticks([min(ESA_V_range):20*beam_energy:max(ESA_V_range)]);
    ylim([-16 16]); yticks([-15:2:15]);
    % ax = gca; ax.YColor = 'w'; ax.XColor = 'w'; h.Color=[1 1 1];
    title('Coincidence Rate')
    set(gca, 'XDir','reverse'); set(gca, 'YDir','normal'); 
    set(gca,'lineWidth',1,'FontSize',12);
    a = gca; b = copyobj(a, gcf);
    set(b,'Xcolor',[1 1 1],'YColor',[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);      
  
    end
    
end

% ===== Parameters for PLOT SAVING ===== 
set(0, 'DefaultFigurePaperPositionMode', 'auto')
set(gcf, 'InvertHardcopy', 'off')
if save_plot == 1
    print('-dpng'  ,'-r300', ['/Users/lykhoo/Library/CloudStorage/OneDrive-PrincetonUniversity/Cal4.0_database/Data Analysis/' ...
        datefld(11:14) '-' datefld(16:17) '-' num2str(dd,'%02d') ...
          '_halfPassband_InnRot_' num2str(inner_angle) '_500eV_Ar_5E-8torr_old_IS_parameters.png'])
else
end
disp('Mission Completed')


