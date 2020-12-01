%% This function sweeps parameter space, running a series of models
%Europa-like conditions
clear;
close all;
addpath core;
seconds_in_year = 3.1558e7;
% Europa
% parameters = struct();
% parameters.Tb = 273;
% parameters.Ts = 100;
% parameters.g  = 1.30;
% parameters.Ro = 1.561e6;
% parameters.Ri = parameters.Ro-4.0e3;     % inner radius of ice shell (m)
% parameters.Rc = parameters.Ro-1.2e5;     % core radius (m)
% parameters.relaxation_parameter = 1e-5;  % europa value
% parameters.tensile_strength = 3e6;
% parameters.perturbation_period = 1e7*seconds_in_year;

% Enceladus
parameters = struct();
parameters.Tb = 273;
parameters.Ts = 100;
parameters.g  = 0.113;
parameters.Ro = 2.52e5;
parameters.Rc = parameters.Ro-1.6e5;     % core radius (m)
parameters.relaxation_parameter = 1e-2;  
parameters.tensile_strength = 3e6;
parameters.perturbation_period = 1e7*seconds_in_year;

ndQ = 6;
dQ = linspace(0.1,0.8,ndQ) ;
nthick = 7;
thicknesses = linspace(4e3,20e3,nthick);
all_results = cell(ndQ,nthick);
all_parameters = cell(ndQ,nthick);
for idQ = 1:length(dQ)
    parfor ithick = 1:length(thicknesses)
        parameters1 = parameters;
        parameters1.Ri = parameters.Ro-thicknesses(ithick);
        parameters1.deltaQonQ = dQ(idQ);
        all_parameters{idQ,ithick} = parameters1;
        all_results{idQ,ithick} = main_cyclic_thermomechanical_model(parameters1);
    end
end
%% Analyze models
all_erupted_volume = zeros(size(all_results));
for ithick=1:nthick
    for idQ=1:ndQ
        results = all_results{idQ,ithick};
        parameters = all_parameters{idQ,ithick};
        ifail = find(results.failure_time>0,1,'last');
        time_mask = results.failure_time >= 1*parameters.perturbation_period/seconds_in_year/1e6;
        % calculate the total volume erupted
        all_erupted_volume(idQ,ithick) = sum( results.failure_erupted_volume(time_mask) );        
    end
end

%% Plotting routines:
for ithick=1:nthick
    for idQ=1:ndQ
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
    contourf(results.time(mask)/seconds_in_year,results.save_depths/1000,results.sigma_t(:,mask),64,'Color','none'); %shading flat;
    hold on
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
        plot(results.failure_time(i)*1e6*[1 1],[results.failure_top(i) results.failure_bottom(i)]/1e3,'r');
    end
    nexttile
    plot(results.time(mask)/seconds_in_year,results.Pex(mask));
    ylabel('Ocean overpressure (Pa)');
    ax2 = gca();
    ax2.Position(3) = ax1.Position(3);
    ax2.XLim = ax1.XLim;
    hold on
    plot(results.failure_time(1:ifail)*1e6,results.failure_P(1:ifail),'ro');
    plot(results.failure_time(1:ifail)*1e6,(results.failure_P(1:ifail)+results.failure_dP(1:ifail)),'g.');
    nexttile
    hold on;
    for i=1:ifail
        plot(results.failure_time(i)*1e6*[1 1],results.failure_erupted_volume(i)/(4*pi*parameters.Ro^2)*[0 1],'b');
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