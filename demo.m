clc;
clear;
close all
font_size = 14;
line_width = 1.2;

%% Parameters
array_num = 20;                 % number of sensors
snapshot_num = 50;              % Snapshots
source_doa = [81.7,95.4,120.7]; % directions of arrival
c = 1520;                       % speed of sound
f = 500;                        % frequency
lambda = c/f;                   % wavelength
d = 0.5*lambda;                 % intersensor spacing
source_num = length(source_doa);% number of signal
snr = -5;                        % SNR
reso = 2;                       % grid resolution
reso_grid = (0:reso:180);       % initial grid

%% Signal generate
rng(1, 'twister' );
S = RandStream('mt19937ar','Seed',5489);
Asignal = exp(-1i*(0:array_num-1)'*2*pi*(d/lambda)*cosd(source_doa));
Xsource = exp(1i*2*pi*rand(source_num,snapshot_num));    % random phase
Ysignal = Asignal*Xsource; 
Y=awgn(Ysignal,snr,'measured');

%% SBL_gh
tic
[Pm_sblgh,~,~] = SBL_gh(Y,d,lambda,reso_grid);
Pm_sblgh = Pm_sblgh/max(Pm_sblgh);
Pm_sblgh = 10*log10(Pm_sblgh);
toc

%% SBL_gaussian
etc = source_num;
tic
[Pm_ogsblgh, search_area]= OGSBL_gh(Y,d,lambda,reso_grid,etc);
Pm_ogsblgh = Pm_ogsblgh/max(Pm_ogsblgh);
Pm_ogsblgh = 10*log10(Pm_ogsblgh);
toc

%% Plot
figure
plot(reso_grid,Pm_sblgh,'LineWidth',line_width);
hold on;
plot(search_area,Pm_ogsblgh,'LineWidth',line_width);
hold on;
plot(source_doa,max(Pm_sblgh),'ro','LineWidth',line_width);
set(gca,'XLim',[0,180],'YLim',[-40,5]);
legend('SBL-GH','OGSBL-GH','True DOAs','fontsize',font_size,'fontname','Times New Roman');
xlabel('\theta (\circ)','fontsize',font_size,'fontname','Times New Roman');
ylabel('Normalized Spectrum(dB)','fontsize',font_size,'fontname','Times New Roman');
hold off;
