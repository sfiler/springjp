% Code for Gain Curve Matrix

%% 00_User Input
close all; clear all; clc;
format long

tdate = '2022.05.12'; % test date
dir = '/Users/lykhoo/Library/CloudStorage/OneDrive-PrincetonUniversity/Documents/Lab/SWAPI Cal/FM-Cal/Data_analysis'
fld = '/HVRAMP-Curves_2022.05.12_H/';

species = 'H';
PCEM_num = [7:1:10]; % PCEM ramp file name for gain curve matrix
SCEM_num = [6:1:9]; % SCEM ramp file name for gain curve matrix

figTitle = 'PCEM-SCEM matrix, ESA = -533V, p+ (8E-6 Torr)';
legendList = {'SCEM +2100V','SCEM +2200V','SCEM +2300V','SCEM +2400V','SCEM +2500V'};
save_folder = '/Analysis';

%% PCEM_SCEM_matrix: alterSCEM
close all; clear data;

figdate = strrep(tdate,'.','_');

for i = 1:length(PCEM_num)
    name = ['HVRAMP_PCEM_',tdate,'_' num2str(PCEM_num(i))];
    path = [dir fld name '.csv'];
    data = importfile_HVramp(path);

    PCEM_V   = table2array(data(:,3)); 
    PCEM_prm = table2array(data(:,8));
    PCEM_sec = table2array(data(:,9));
    PCEM_coi = table2array(data(:,10));

    % keep the sort index in "sortIdx"
    [PCEM_V,sortIdx] = sort(PCEM_V,'descend');
    % sort B using the sorting index
    PCEM_prm_bin(:,i) = PCEM_prm(sortIdx);
    PCEM_sec_bin(:,i) = PCEM_sec(sortIdx);
    PCEM_coi_bin(:,i) = PCEM_coi(sortIdx);
    PCEM_prm_plot(:,i) = PCEM_prm_bin(:,i);
    PCEM_sec_plot(:,i) = PCEM_sec_bin(:,i);
    PCEM_coi_plot(:,i) = PCEM_coi_bin(:,i);
    epsilon_D_bin(:,i) = (PCEM_coi_plot(:,i).*PCEM_coi_plot(:,i))./(PCEM_prm_plot(:,i).*PCEM_sec_plot(:,i));
    PCEM_Vplot = PCEM_V;
    figure(1)
    h1 = plot(PCEM_Vplot,epsilon_D_bin(:,i),'-o','linewidth',2.0); hold on;
    set(h1, 'MarkerFaceColor', get(h1,'Color'));
    xlim([-2400 -1700]); xticks([-2400:50:-1700])
    ylim([0   .1]); yticks([0:.01:.1])
    xlabel('Primary CEM voltage (V)')
    ylabel(['Absolute detection efficiency (' char(949) '_D)'])
    legend(legendList,'location','southeast')
    title(figTitle);
    set(gca, 'XDir','reverse'); grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder '/HVGain_matrix_PCEMSweep_' species '_ADE_',figdate,'.png'])
    end

    figure(2)
    h2 = plot(PCEM_Vplot,PCEM_prm_plot(:,i)./PCEM_sec_plot(:,i),'-o','linewidth',2.0); hold on;
    set(h2, 'MarkerFaceColor', get(h2,'Color'));
    xlim([-2400 -1700]); xticks([-2400:50:-1700])
    ylim([0   1.4]);   yticks([0:.1:1.4])
    xlabel('Primary CEM voltage (V)')
    ylabel(['Primary/Secondary Ratio'])
    legend(legendList,'location','southeast')
    title(figTitle);
    set(gca, 'XDir','reverse'); grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder '/HVGain__13_matrix_PCEMSweep_' species '_PSR_',figdate,'.png'])
    end

    figure(3)
    h3 = plot(PCEM_Vplot,PCEM_prm_plot(:,i),'-o','linewidth',2.0); hold on;
    set(h3, 'MarkerFaceColor', get(h3,'Color'));
    xlim([-2400 -1700]); xticks([-2400:50:-1700])
%     ylim([0   500]);   yticks([0.2:.1:1.2])
    xlabel('Primary CEM voltage (V)')
    ylabel(['PCEM Count Rate (Hz)'])
    legend(legendList,'location','southeast')
    title(figTitle);
    set(gca, 'XDir','reverse'); grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i ==4
    print('-dpng'  ,'-r300', [dir save_folder '/HVGain__13_matrix_PCEMSweep_' species '_Pcounts_',figdate,'.png'])
    end

    figure(4)
    h4 = plot(PCEM_Vplot,PCEM_sec_plot(:,i),'-o','linewidth',2.0); hold on;
    set(h4, 'MarkerFaceColor', get(h4,'Color'));
    xlim([-2400 -1700]); xticks([-2400:50:-1700])
    % ylim([0.2   1.2]);   yticks([0.2:.1:1.2])
    xlabel('Primary CEM voltage (V)')
    ylabel(['SCEM Count Rate (Hz)'])
    legend(legendList,'location','southeast')
    title(figTitle);
    set(gca, 'XDir','reverse'); grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder '/HVGain__13_matrix_PCEMSweep_' species '_Scounts_',figdate,'.png'])
    end

    figure(5)
    h5 = plot(PCEM_Vplot,PCEM_coi_plot(:,i),'-o','linewidth',2.0); hold on;
    set(h5, 'MarkerFaceColor', get(h5,'Color'));
    xlim([-2400 -1700]); xticks([-2400:50:-1700])
    % ylim([0.2   1.2]);   yticks([0.2:.1:1.2])
    xlabel('Primary CEM voltage (V)')
    ylabel(['Coincidence Rate (Hz)'])
    legend(legendList,'location','southeast')
    title(figTitle)
    set(gca, 'XDir','reverse'); grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder  '/HVGain__13_matrix_PCEMSweep_' species '_Ccounts_',date,'.png'])
    end


end
disp('Mission Completed')



%% PCEM_SCEM_matrix: alterPCEM
close all; clear data epsilon_D_bin;

for i = 1:length(SCEM_num)
    name = ['HVRAMP_SCEM_',tdate,'_' num2str(SCEM_num(i))];
    path = [dir fld name '.csv'];
    data = importfile_HVramp(path);
    
    SCEM_V   = table2array(data(:,3)); 
    SCEM_prm = table2array(data(:,8));
    SCEM_sec = table2array(data(:,9));
    SCEM_coi = table2array(data(:,10));
    % keep the sort index in "sortIdx"
    [SCEM_V,sortIdx] = sort(SCEM_V);
    % sort B using the sorting index
    SCEM_prm_bin = SCEM_prm(sortIdx);
    SCEM_sec_bin = SCEM_sec(sortIdx);
    SCEM_coi_bin = SCEM_coi(sortIdx);
    SCEM_prm_plot = SCEM_prm_bin;
    SCEM_sec_plot = SCEM_sec_bin;
    SCEM_coi_plot = SCEM_coi_bin;
    SCEM_V_plot = SCEM_V;
    epsilon_D_bin = (SCEM_coi_plot.*SCEM_coi_plot)./(SCEM_prm_plot.*SCEM_sec_plot);
    
    figure(1)
    h1 = plot(SCEM_V_plot,epsilon_D_bin,'-o','linewidth',2.0); hold on;
    set(h1, 'MarkerFaceColor', get(h1,'Color'));
    xlim([1700 2500]); xticks([1700:50:2500])
    ylim([0   .1]); yticks([0:.01:.1])
    xlabel('Secondary CEM voltage (V)')
    ylabel(['Absolute detection efficiency (' char(949) '_D)'])
    legend(legendList,'location','southeast'); 
    title(figTitle);
    % set(gca, 'XDir','reverse'); 
    grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder  '/HVGain__13_matrix_SCEMSweep_' species '_ADE_',figdate,'.png'])
    end
    
    figure(2)
    h2 = plot(SCEM_V_plot,SCEM_sec_plot./SCEM_prm_plot,'-o','linewidth',2.0); hold on;
    set(h2, 'MarkerFaceColor', get(h2,'Color'));
    xlim([1700 2500]); xticks([1700:50:2500])
    ylim([0   1.6]);   yticks([0:.1:1.6])
    xlabel('Secondary CEM voltage (V)')
    ylabel(['Secondary/Primary Ratio'])
    legend(legendList,'location','southeast'); 
    title(figTitle);
    % set(gca, 'XDir','reverse'); 
    grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder  '/HVGain__13_matrix_SCEMSweep_' species '_SPR_',figdate,'.png'])
    end
    
    figure(3)
    h3 = plot(SCEM_V_plot,SCEM_prm_plot,'-o','linewidth',2.0); hold on;
    set(h3, 'MarkerFaceColor', get(h3,'Color'));
    xlim([1700 2500]); xticks([1700:50:2500])
    % ylim([0   1.4]);   yticks([0:.1:1.4])
    xlabel('Secondary CEM voltage (V)')
    ylabel(['PCEM Count Rate (Hz)'])
    legend(legendList,'location','southeast'); 
    title(figTitle);
    % set(gca, 'XDir','reverse'); 
    grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder  '/HVGain__13_matrix_SCEMSweep_' species '_Pcounts_',figdate,'.png'])
    end
    
    figure(4)
    h4 = plot(SCEM_V_plot,SCEM_sec_plot,'-o','linewidth',2.0); hold on;
    set(h4, 'MarkerFaceColor', get(h4,'Color'));
    xlim([1700 2500]); xticks([1700:50:2500])
    % ylim([0   1.4]);   yticks([0:.1:1.4])
    xlabel('Secondary CEM voltage (V)')
    ylabel(['SCEM Count Rate (Hz)'])
    legend(legendList,'location','southeast'); 
    title(figTitle);
    % set(gca, 'XDir','reverse'); 
    grid on;
    a = gca; a.XTickLabelRotation = 75;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder  '/HVGain__13_matrix_SCEMSweep_' species '_Scounts_',figdate,'.png'])
    end
    
    figure(5)
    h5 = plot(SCEM_V_plot,SCEM_coi_plot,'-o','linewidth',2.0); hold on;
    set(h5, 'MarkerFaceColor', get(h5,'Color'));
    xlim([1700 2500]); xticks([1700:50:2500])
    % ylim([0   1.4]);   yticks([0:.1:1.4])
    xlabel('Secondary CEM voltage (V)')
    ylabel(['Coincidence Rate (Hz)'])
    legend(legendList,'location','southeast'); 
    title(figTitle);
    % set(gca, 'XDir','reverse'); 
    grid on;
    set(gca,'lineWidth',1,'FontSize',14); %,'FontWeight','bold'
    a = gca; a.XTickLabelRotation = 75;
    if i == 4
    print('-dpng'  ,'-r300', [dir save_folder '/HVGain__13_matrix_SCEMSweep_' species '_Ccounts_',figdate,'.png'])
    end

end

disp('Mission Completed')



