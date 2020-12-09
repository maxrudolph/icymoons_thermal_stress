%% This function sweeps parameter space, running a series of models
%Europa-like conditions
clear;
close all;
addpath core;
seconds_in_year = 3.1558e7;

for moon=0:1
    if moon==0
        % Europa
        parameters = struct();
        parameters.Tb = 273;
        parameters.Ts = 100;
        parameters.g  = 1.30;
        parameters.Ro = 1.561e6;
        parameters.Rc = parameters.Ro-1.2e5;     % core radius (m)
        parameters.relaxation_parameter = 1e-4;  % europa value
        parameters.tensile_strength = 3e6;
        parameters.perturbation_period = 1e8*seconds_in_year;
        parameters.save_start = parameters.perturbation_period*5;
        parameters.save_interval = parameters.perturbation_period/1000;
        parameters.end_time = parameters.perturbation_period*10;
        parameters.label = 'Europa';
    else
        % Enceladus
        parameters = struct();
        parameters.Tb = 273;
        parameters.Ts = 100;
        parameters.g  = 0.113;
        parameters.Ro = 2.52e5;
        parameters.Rc = parameters.Ro-1.6e5;     % core radius (m)
        parameters.relaxation_parameter = 1e-2;
        parameters.tensile_strength = 3e6;
        parameters.perturbation_period = 1e8*seconds_in_year;
        parameters.save_start = parameters.perturbation_period*5;
        parameters.save_interval = parameters.perturbation_period/1000;
        parameters.end_time = parameters.perturbation_period*10;
        parameters.label = 'Enceladus';
    end
    
    ndQ = 15;
    dQ = linspace(0.1,0.8,ndQ) ;
    nthick = 33;
    thicknesses = linspace(4e3,20e3,nthick);
    all_results = cell(ndQ,nthick);
    all_parameters = cell(ndQ,nthick);
    
    for idQ = 1:length(dQ)
        for ithick = 1:length(thicknesses)
            parameters1 = parameters;
            parameters1.Ri = parameters.Ro-thicknesses(ithick);
            parameters1.deltaQonQ = dQ(idQ);
            all_parameters{idQ,ithick} = parameters1;
        end
    end
    parfor irun=1:ndQ*nthick
        all_results{irun} = main_cyclic_thermomechanical_model(all_parameters{irun});
    end
    save([parameters.label '_workspace.mat'],'all_parameters','all_results','ndQ','nthick','thicknesses','dQ','-v7.3');
end
%% Load results
clearvars -except seconds_in_year;
seconds_in_year = 3.1558e7;

load('Enceladus_workspace.mat');
% load('Europa_workspace.mat');

%% Analyze models
all_erupted_volume = zeros(size(all_results));
all_failure_events = zeros(size(all_results));
for ithick=1:nthick
    for idQ=1:ndQ
        results = all_results{idQ,ithick};
        parameters = all_parameters{idQ,ithick};
        ifail = find(results.failure_time>0,1,'last');
        time_mask = results.failure_time >= 2*parameters.perturbation_period/seconds_in_year/1e6;
        % calculate the total volume erupted
        all_erupted_volume(idQ,ithick) = sum( results.failure_erupted_volume(time_mask) );
        all_failure_events(idQ,ithick) = nnz( results.failure_time(time_mask) );
    end
end
figure();
addpath ~/sw/matlab/crameri

t=tiledlayout(2,1,'Padding','none','TileSpacing','none','Units','centimeters');
t.OuterPosition(3:4) = [9 9];
nexttile;
dt = max(results.time) - 2*parameters.perturbation_period;
contourf(thicknesses/1e3,dQ,all_erupted_volume/(4*pi*parameters.Ro^2)/(dt/seconds_in_year/1e6),16,'Color','none');
hcb=colorbar();
hcb.Label.String = 'Erupted Volume (m/Myr)';
colormap(crameri('davos'));
ylabel('\Delta q/q_0');
xlabel('Average Thickness (km)');
if( max(all_erupted_volume(:))==0 )
    caxis([0 1])
end
title(parameters.label);
nexttile;
dt = max(results.time) - 2*parameters.perturbation_period;
contourf(thicknesses/1e3,dQ,all_failure_events/(dt/parameters.perturbation_period),16,'Color','none');
hcb=colorbar();

hcb.Label.String = 'Failure events per cycle';
ylabel('\Delta q/q_0');
xlabel('Average Thickness (km)');
colormap(crameri('davos'));
exportgraphics(gcf,[parameters.label '_regimes.eps']);

%% Look at phase lag between qb and thickness
all_phase_lags = zeros(ndQ,nthick);
all_ampfrac = zeros(ndQ,nthick);
for ithick=1:nthick
    for idQ=1:ndQ
        
        rho_i = 900;
        Cp = 2100;              % heat capacity of ice (J/kg/K)
        kappa = 1e-6;           % Thermal diffusivity (m/s/s)
        k=kappa*rho_i*Cp;
        
        results = all_results{idQ,ithick};
        parameters = all_parameters{idQ,ithick};       
        
        actual_dz = (results.Ri-results.z) - mean(results.Ri-results.z);
        actual_thickness = parameters.Ro - (results.Ri-results.z);
        %         plot(results.time/seconds_in_year, (actual_dz/std(actual_dz)) );
%         figure, plot(results.time/seconds_in_year, actual_thickness ); hold on
        ss_thickness = k*(273-100)./results.qb;
        [ss_peaks,ss_ind] = findpeaks(ss_thickness);
        [actual_peaks,actual_ind] = findpeaks(actual_thickness);
        phase_lag = median(results.time(actual_ind) - results.time(ss_ind))/seconds_in_year/1e6;
        all_phase_lags(idQ,ithick) = phase_lag;
        all_ampfrac(idQ,ithick) = (max(actual_thickness)-min(actual_thickness))/(max(ss_thickness)-min(ss_thickness));
%         plot(results.time/seconds_in_year,ss_thickness,'--');
    end
end
figure;
subplot(2,1,1);
contourf(thicknesses/1e3,dQ,all_phase_lags,16,'Color','none'); colorbar();
ylabel('\Delta q/q_0');
xlabel('Average Thickness (km)');
subplot(2,1,2);
contourf(thicknesses/1e3,dQ,all_ampfrac,16,'Color','none'); colorbar();
ylabel('\Delta q/q_0');
xlabel('Average Thickness (km)');
%% Plotting routines:
for ithick=1:1
    for idQ=1:3:ndQ
        results = all_results{idQ,ithick};
        parameters = all_parameters{idQ,ithick};
        isave = find(results.time>0,1,'last');
        ifail = find(results.failure_time>0,1,'last');
        mask = 1:(isave);
        
        figure();
        subplot(3,1,1);
        plot(results.time(mask)/seconds_in_year,results.z(mask))
        ylabel('Amount of freezing (m)');
        subplot(3,1,2);
        plot(results.time(mask)/seconds_in_year,results.sigma_t(1,mask),'DisplayName',sprintf('%.02f km',results.save_depths(1)/1000));
        hold on
        plot(results.time(mask)/seconds_in_year,results.sigma_t(10,mask),'DisplayName',sprintf('%.02f km',results.save_depths(10)/1000));
        plot(results.time(mask)/seconds_in_year,results.sigma_t(20,mask),'DisplayName',sprintf('%.02f km',results.save_depths(20)/1000));
        legend();
        ylabel('\sigma_t (Pa)')
        subplot(3,1,3);
        plot(results.time(mask)/seconds_in_year,results.Pex(mask));
        ylabel('Overpressure (Pa)');
        
        %% Pseudocolor stress plot
        figure();
        t=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
        nexttile
        contourf(results.time(mask)/seconds_in_year,results.save_depths/1000,results.sigma_t(:,mask),16,'Color','none'); %shading flat;
        hold on
        caxis([-1 1]*max(abs(caxis())));
        colormap(crameri('vik'));
        plot(results.time(mask)/seconds_in_year,((parameters.Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        set(gca,'YLim',[0 ceil(1+max(((parameters.Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        ax1 = gca();
        hcb = colorbar();
        hcb.Label.String = 'Tensile Stress (Pa)';
        xlabel('Time (years)');
        ylabel('Depth (km)');
        hold on;
        for i=1:ifail
            plot(results.failure_time(i)/seconds_in_year*[1 1],[results.failure_top(i) results.failure_bottom(i)]/1e3,'r');
        end
        title(sprintf('\\Delta Q/Q_0 = %.2f, H_0=%.1f km',parameters.deltaQonQ,(parameters.Ro-parameters.Ri)/1000));
        
        nexttile
        plot(results.time(mask)/seconds_in_year,results.Pex(mask));
        ylabel('Ocean overpressure (Pa)');
        ax2 = gca();
        ax2.Position(3) = ax1.Position(3);
        ax2.XLim = ax1.XLim;
        hold on
        plot(results.failure_time(1:ifail)/seconds_in_year,results.failure_P(1:ifail),'ro');
        plot(results.failure_time(1:ifail)/seconds_in_year,(results.failure_P(1:ifail)+results.failure_dP(1:ifail)),'g.');
        nexttile
        hold on;
        for i=1:ifail
            plot(results.failure_time(i)/seconds_in_year*[1 1],results.failure_erupted_volume(i)/(4*pi*parameters.Ro^2)*[0 1],'b');
            %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_volumechange(i)/(4*pi*Ro^2),'go');
            %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_pressurechange(i)/(4*pi*Ro^2),'rx');
        end
        ylabel('Erupted volume (m)');
        xlabel('Time (years)');
        ax3=gca();
        ax3.XLim = ax1.XLim;
        ax3.Position(3) = ax1.Position(3);
        ax3.Box = 'on';
        fig = gcf();
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition(3) = 12.00;
        fig.Color = 'w';
        %     exportgraphics(gcf,'test.eps','ContentType','vector');
    end
end
% End main code
