%% This function sweeps parameter space, running a series of models
%Europa-like conditions
clear;
close all;
addpath core;
seconds_in_year = 3.1558e7;
do_runs = true
strength_label = '1MPa';

if do_runs
    for moon=1:1
        if moon==0
            % Europa
            parameters = struct();
            parameters.g  = 1.30;
            parameters.Ro = 1.561e6;
            parameters.Rc = parameters.Ro-1.2e5;     % core radius (m)
            parameters.relaxation_parameter = 1e-4;  % europa value            
            parameters.label = 'Europa';
        else
            % Enceladus
            parameters = struct();            
            parameters.g  = 0.113;
            parameters.Ro = 2.52e5;
            parameters.Rc = parameters.Ro-1.6e5;     % core radius (m)
            parameters.relaxation_parameter = 1e-1;           
            parameters.label = 'Enceladus';
        end
        
        parameters.tensile_strength = 1e6;
        parameters.viscosity_model = 0;  % 0 = Nimmo (2004), 1 = Goldsby and Kohlstedt (2001)
        parameters.nr = 512;
        parameters.Tb = 273;
        parameters.Ts = 100;
        parameters.k = @(T) 651./T; % Petrenko and Whitworth 1999
        parameters.perturbation_period = 1e8*seconds_in_year;
        parameters.save_start = parameters.perturbation_period*5;
        parameters.save_interval = parameters.perturbation_period/1000;
        parameters.end_time = parameters.perturbation_period*10;
        

          ndQ = 4;
          dQ = linspace(0.1,0.8,ndQ);
          nthick = 3;
          thicknesses = logspace(log10(2e3),log10(20e3),nthick);
%         ndQ = 1;
%         nthick = 1;
%         dQ = [0.8];
%         thicknesses = [2e4];
        all_results = cell(ndQ,nthick);
        all_parameters = cell(ndQ,nthick);
        
        for idQ = 1:length(dQ)
            for ithick = 1:length(thicknesses)
                parameters1 = parameters;
                parameters1.Ri = parameters.Ro-thicknesses(ithick);
                parameters1.deltaQonQ = dQ(idQ);
                if thicknesses(ithick) < 2e3
                    parameters1.relaxation_parameter = parameters1.relaxation_parameter/10;
                end
                parameters1.nr = min([512 ceil(thicknesses(ithick)/20)]);
                all_parameters{idQ,ithick} = parameters1;
            end
        end
        % parfor here for real runs:
        for irun=4:ndQ*nthick
            all_results{irun} = main_cyclic_thermomechanical_model(all_parameters{irun});
        end
        save([parameters.label '_' num2str(parameters.tensile_strength/1e6) 'MPa_workspace.mat'],'all_parameters','all_results','ndQ','nthick','thicknesses','dQ','-v7.3');
    end
else% postprocess:
    %% Load results
    close all;
    addpath ~/sw/matlab/crameri
    
    for moon=0:1
        clearvars -except seconds_in_year moon strength_label;
        seconds_in_year = 3.1558e7;
        disp('loading data');
        if moon==0
            load(['./Europa_' strength_label '_workspace.mat']);
        else
            load(['./Enceladus_' strength_label '_workspace.mat']);
        end
        
        % data selection
        % to make plotting simpler, select data corresponding to
        % interesting part of parameter space for each tensile strength
        %         if strcmp(strength_label,'1MPa')
        %             mask = thicknesses <= 1.1e4 & thicknesses >= 3e3;
        %         else
        %             mask = thicknesses <= 2e4   & thicknesses >= 3e3;
        %         end
        %         all_results = all_results(:,mask);
        %         all_parameters = all_parameters(:,mask);
        %         thicknesses = thicknesses(mask);
        %         nthick = length(thicknesses);
        
        %% Analyze models
        disp('analyzing models');
        all_erupted_volume = zeros(size(all_results));
        all_failure_events = zeros(size(all_results));
        all_failure_fraction = zeros(size(all_results));
        all_actual_thickness = zeros(size(all_results));
        all_failure_initial = zeros(size(all_results));
        all_first_eruption_time = zeros(ndQ,nthick);
        all_reach_ocean = zeros(ndQ,nthick);
        
        % Convert thickness into Q0
        Qtots = zeros(size(thicknesses));
        
        for ithick=1:nthick
            for idQ=1:ndQ
                results = all_results{idQ,ithick};
                parameters = all_parameters{idQ,ithick};
                if idQ == 1
                    [~,T_last,qr] = find_steady_T(parameters.Ri,parameters.Ro,parameters.Tb,parameters.Ts,[parameters.Ri parameters.Ro]);
                    Q0 = 4*pi*parameters.Ro^2 * qr(end);
                    Qtots(ithick)=Q0;
                end
                
                
                ifail = find(results.failure_time>0,1,'last');
                results.Qtot = results.Qtot(1:length(results.qb));
                all_results{idQ,ithick} = results;
                time_mask = results.failure_time >= 2*parameters.perturbation_period/seconds_in_year/1e6;
                % calculate the total volume erupted
                if length(results.failure_erupted_volume) < length(time_mask)
                    warning('padding failure_erupted_volume');
                    results.failure_erupted_volume(end+1:length(time_mask)) = 0;
                    all_results{idQ,ithick} = results;
                end
                if ~isempty(time_mask)
                    all_erupted_volume(idQ,ithick) = sum( results.failure_erupted_volume(time_mask) );
                    all_failure_events(idQ,ithick) = nnz( results.failure_time(time_mask) );
                    all_failure_initial(idQ,ithick) = median( results.failure_initial(time_mask) );
                else
                    all_erupted_volume(idQ,ithick) = NaN;
                    all_failure_events(idQ,ithick) = NaN;
                    all_failure_initial(idQ,ithick) = NaN;
                end
                
                
                % calculate the average thickness.
                thickness = parameters.Ro-(results.Ri-results.z);
                average_thickness = 1/(results.time(end)-results.time(1))*sum( (thickness(1:end-1) + thickness(2:end))/2 .* diff(results.time));
                all_actual_thickness(idQ,ithick) = average_thickness;
                % approximate ice shell thickness at time of failure
                if any(results.failure_thickness)
                    failure_thickness = interp1(results.time,parameters.Ro-(results.Ri-results.z),results.failure_time(time_mask));
                    all_failure_fraction(idQ,ithick) = max(results.failure_thickness ./ failure_thickness);
                    all_reach_ocean(idQ,ithick) = max(results.failure_bottom ./failure_thickness);
                else
                    all_failure_fraction(idQ,ithick) = NaN;
                    all_reach_ocean(idQ,ithick) = NaN;
                end
                % if there are eruptions, calculate timing relative to
                % thickening
                
                if ~isempty(results.failure_eruption_time)
                    % for each eruption time, find preceding ice shell thickness minimum.
                    [~,p1] = findpeaks(-thickness,'MinPeakProminence',100); % find thickness minima (anti-peaks)
                    thickness_minima = results.time(p1);
                    thickness_minima_mod = mod(thickness_minima,parameters.perturbation_period);
                    eruption_time_mod = mod(results.failure_eruption_time,parameters.perturbation_period);
                    first_eruption_time = min(eruption_time_mod) - mean(thickness_minima_mod);
                    if first_eruption_time <0
                        first_eruption_time= first_eruption_time+parameters.perturbation_period;
                    end
                    all_first_eruption_time(idQ,ithick) = first_eruption_time/seconds_in_year/1e6;
                    
                else
                    all_first_eruption_time(idQ,ithick) = NaN;
                end
            end
        end
        %%
        disp('plotting');
        if moon==0
            figure(901);
            clf;
            t1=tiledlayout(3,2,'Padding','none','TileSpacing','compact','Units','centimeters');
            t1.OuterPosition(3:4) = [18 12];
        else
            figure(901);
        end
        
        
        %         nexttile(1+moon);
        %         dt = max(results.time) - 2*parameters.perturbation_period;
        %         contourf(thicknesses/1e3,dQ,all_erupted_volume/(4*pi*parameters.Ro^2)/(dt/seconds_in_year/1e6),16,'Color','none');
        %         hcb=colorbar();
        %         hcb.Label.String = 'Erupted Volume (m/Myr)';
        %         colormap(crameri('davos'));
        %         if ~moon
        %             ylabel('\Delta q/q_0');
        %         end
        %         %         xlabel('Equilibrium Thickness (km)');
        %         if( max(all_erupted_volume(:))==0 )
        %             caxis([0 1])
        %         end
        %         title(parameters.label);
        %         text(0.05,0.85,char(65+4*moon),'FontSize',14,'Units','normalized');
        
        % Panel 2
        nexttile(1+moon);
        contourf(thicknesses/1e3,dQ,all_failure_fraction,20,'Color','none');
        hold on
        contour(thicknesses/1e3,dQ,all_failure_fraction,0.985*[1 1],'k');
        contour(thicknesses/1e3,dQ,all_reach_ocean,[1 1]*0.985,'r--');
        caxis([0 1]);
        hcb = colorbar();
        hcb.Label.String = 'Fractional penetration';
        if ~moon
            ylabel('\Delta q/q_0');
        end
        text(0.05,0.85,char(65+3*moon+0),'FontSize',14,'Units','normalized');
        set(gca,'Color',[1 1 1]*0.85)
        
        %         xlabel('Equilibrium Thickness (km)');
        
        % Panel 3
        nexttile(3+moon);
        contourf(thicknesses/1e3,dQ,all_failure_initial/1e3,20,'Color','none');
        hcb=colorbar();
        hcb.Label.String = 'Failure Depth (km)';
        if ~moon
            ylabel('\Delta q/q_0');
        end
        text(0.05,0.85,char(65+3*moon+1),'FontSize',14,'Units','normalized');
        set(gca,'Color',[1 1 1]*0.85)
        % Panel 4
        nexttile(5+moon);
        dt = max(results.time) - 2*parameters.perturbation_period;
        failures_per_cycle = all_failure_events/(dt/parameters.perturbation_period);
        contourf(thicknesses/1e3,dQ,failures_per_cycle,linspace(0,200,20),'Color','none');
        hcb=colorbar();
        text(0.05,0.85,char(65+3*moon+2),'FontSize',14,'Units','normalized');
        set(gca,'Color',[1 1 1]*0.85)
        hcb.Label.String = 'Cracks per cycle';
        if ~moon
            ylabel('\Delta q/q_0');
        end
        
        xlabel('Equilibrium Thickness (km)');
        colormap(crameri('davos'));
        
        if moon == 1
            exportgraphics(gcf,['Figure2_' strength_label '_regimes.eps']);
        end
        %% Look at phase lag between qb and thickness
        all_phase_lags = zeros(ndQ,nthick);
        all_ampfrac = zeros(ndQ,nthick);
        all_Qbar = zeros(ndQ,nthick);
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
                %                 ss_thickness = k*(273-100)./results.qb;
                
                qsurf = results.Qtot./(4*pi*parameters.Ro^2);
                ss_thickness = find_steady_h(parameters.Ro,parameters.Tb,parameters.Ts,qsurf);
                if any( imag(ss_thickness) ~= 0 )
                    ss_thickness = real(ss_thickness);
                    warning('complex thickness encountered - this case is not reasonable');
                end
                % old method for estimating phase lag:
                % [ss_peaks,ss_ind] = findpeaks(ss_thickness,'MinPeakProminence',std(ss_thickness)/4);
                % [actual_peaks,actual_ind] = findpeaks(actual_thickness,'MinPeakProminence',std(actual_thickness)/4);
                % phase_lag = median(results.time(actual_ind) - results.time(ss_ind))/seconds_in_year/1e6;
                % interpolate onto uniform time vector
                % Estimate the phase lag using MATLAB builtin finddelay:
                t_tmp = linspace(results.time(1),results.time(end),length(results.time)); % uniform time vector
                dt = t_tmp(2) - t_tmp(1);
                actual_thickness_u = interp1(results.time,actual_thickness,t_tmp); % actual thickness, uniform time vector.
                ss_thickness_u = interp1(results.time,ss_thickness,t_tmp); % steady state thickness, uniform time vector
                %                 phase_lag = finddelay(ss_thickness_u,actual_thickness_u)*dt/seconds_in_year/1e6;
                [~,p1] = findpeaks(actual_thickness_u,'MinPeakProminence',100);
                [~,p2] = findpeaks(ss_thickness_u,'MinPeakProminence',100);
                phase_lag = mean( abs(t_tmp(p2) - t_tmp(p1)) )/seconds_in_year/1e6;
                
                all_phase_lags(idQ,ithick) = phase_lag;
                all_ampfrac(idQ,ithick) = (max(actual_thickness)-min(actual_thickness))/(max(ss_thickness)-min(ss_thickness));
                % calculate average heat flow
                tmp = cumtrapz(results.time,results.Qtot)/(results.time(end)-results.time(1));
                all_Qbar(idQ,ithick) = tmp(end);
                %         plot(results.time/seconds_in_year,ss_thickness,'--');
            end
        end
        if moon==0
            figure(902);
            clf;
            t1=tiledlayout(3,2,'Padding','none','TileSpacing','compact','Units','centimeters');
            t1.OuterPosition(3:4) = [18 12];
        else
            figure(902);
        end
        
        nexttile(1+moon)
        contourf(thicknesses/1e3,dQ,all_phase_lags,20,'Color','none'); crameri('davos')
        hcb=colorbar();
        hcb.Label.String='Lag (Myr)';
        if ~moon
            ylabel('\Delta q/q_0');
        end
        %         xlabel('Equilibrium Thickness (km)');
        text(0.05,0.85,char('A'+3*moon),'Units','normalized','FontSize',14);
        title(parameters.label);
        
        nexttile(3+moon)
        contourf(thicknesses/1e3,dQ,all_ampfrac,20,'Color','none'); crameri('davos')
        hcb=colorbar();
        hcb.Label.String='\Delta h/\Delta h_{eq}';
        if ~moon
            ylabel('\Delta q/q_0');
        end
        %         xlabel('Equilibrium Thickness (km)');
        text(0.05,0.85,char('B'+3*moon),'Units','normalized','FontSize',14);
        caxis([0 1])
        
        nexttile(5+moon)
        contourf(thicknesses/1e3,dQ,all_actual_thickness/1e3,20,'Color','none'); crameri('davos')
        hcb=colorbar();
        hcb.Label.String='Average Thickness (km)';
        if ~moon
            ylabel('\Delta q/q_0');
        end
        xlabel('Equilibrium Thickness (km)');
        text(0.05,0.85,char('C'+3*moon),'Units','normalized','FontSize',14);
        if moon == 1
            exportgraphics(gcf,['Figure3_combined_lag_eqlm_' strength_label '.eps']);
        end
        %% plot total heat flow
        if moon==0
            figure(905);
            clf;
            t2 = tiledlayout(1,2,'Padding','none','TileSpacing','none','Units','centimeters');
            t2.OuterPosition(3:4) = [18 4];
        else
            figure(905);
        end
        nexttile(1+moon);
        contourf(thicknesses/1e3,dQ,all_Qbar/1e12,20,'Color','none'); crameri('davos')
        hcb=colorbar();
        hcb.Label.String='Total Heat Flow (TW)';
        title(parameters.label);
        ylabel('\Delta q/q_0');
        xlabel('Equilibrium Thickness (km)');
        if moon == 1
            exportgraphics(gcf,['Figure4_' strength_label '.eps']);
        end
        %% Plot maximum fractional penetration
        figure();
        
        
        %% Plot outcomes of individual runs
        run_plots = true;
        if run_plots
            for ithick = [1 nthick]
                for idQ = [1 ndQ]
                    %             for ithick=[1 3 14 20 nthick]
                    %                 for idQ=[ 5 ndQ]
                    results = all_results{idQ,ithick};
                    parameters = all_parameters{idQ,ithick};
                    isave = find(results.time>0,1,'last');
                    ifail = find(results.failure_time>0,1,'last');
                    mask = 1:(isave);
                    
                    %                 figure();
                    %                 subplot(3,1,1);
                    %                 plot(results.time(mask)/seconds_in_year,results.z(mask))
                    %                 ylabel('Amount of freezing (m)');
                    %                 subplot(3,1,2);
                    %                 plot(results.time(mask)/seconds_in_year,results.sigma_t(1,mask),'DisplayName',sprintf('%.02f km',results.save_depths(1)/1000));
                    %                 hold on
                    %                 plot(results.time(mask)/seconds_in_year,results.sigma_t(10,mask),'DisplayName',sprintf('%.02f km',results.save_depths(10)/1000));
                    %                 plot(results.time(mask)/seconds_in_year,results.sigma_t(20,mask),'DisplayName',sprintf('%.02f km',results.save_depths(20)/1000));
                    %                 legend();
                    %                 ylabel('\sigma_t (Pa)')
                    %                 subplot(3,1,3);
                    %                 plot(results.time(mask)/seconds_in_year,results.Pex(mask));
                    %                 ylabel('Overpressure (Pa)');
                    
                    %% Pseudocolor stress plot
                    figure();
                    t=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
                    nexttile
                    
                    sigplot = results.sigma_t(:,mask);
                    for i=1:size(results.sigma_t,2)
                        sigplot(:,i) = sigplot(:,i) - 0*rho_i*parameters.g*results.save_depths'; % optionally subtract hydrostatid pressure (compression negative)
                    end
                    contourf(results.time(mask)/seconds_in_year/1e6,results.save_depths/1000,sigplot,16,'Color','none'); %shading flat;
                    hold on
                    caxis([-1 1]*max(abs(caxis())));
                    colormap(crameri('vik'));
                    plot(results.time(mask)/seconds_in_year/1e6,((parameters.Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
                    % plot the location of initial failure.
                    plot(results.failure_time/seconds_in_year/1e6,results.failure_initial/1000,'r.')
                    set(gca,'YLim',[0 ceil(1+max(((parameters.Ro-results.Ri(mask))+results.z(mask))/1000))]);
                    set(gca,'YDir','reverse');
                    ax1 = gca();
                    hcb = colorbar();
                    hcb.Label.String = 'Tensile Stress (Pa)';
                    %                     xlabel('Time (Myr)');
                    ylabel('Depth (km)');
                    hold on;
                    for i=1:ifail
                        plot(results.failure_time(i)/1e6/seconds_in_year*[1 1],[results.failure_top(i) results.failure_bottom(i)]/1e3,'r');
                    end
                    title(sprintf('%s : \\Delta Q/Q_0 = %.2f, H_0=%.1f km',parameters.label,parameters.deltaQonQ,(parameters.Ro-parameters.Ri)/1000));
                    
                    % Plot pressure time series
                    nexttile
                    plot(results.time(mask)/seconds_in_year/1e6,results.Pex(mask));
                    ylabel('Ocean overpressure (Pa)');
                    ax2 = gca();
                    ax2.Position(3) = ax1.Position(3);
                    ax2.XLim = ax1.XLim;
                    hold on
                    plot(results.failure_time(1:ifail)/seconds_in_year/1e6,results.failure_P(1:ifail),'ro');
                    plot(results.failure_time(1:ifail)/seconds_in_year/1e6,(results.failure_P(1:ifail)+results.failure_dP(1:ifail)),'g.');
                    thickness = parameters.Ro - (results.Ri-results.z);
                    b = thickness*(1-rho_i/1000); % depth of neutral buoyancy
                    
                    %                     plot(results.time(mask)/seconds_in_year/1e6,1000*parameters.g*b,'k');
                    plot(results.time(mask)/seconds_in_year/1e6,results.Pex_crit(mask),'k','LineWidth',1);
                    nexttile
                    hold on;
                    % compute whether the pressure at the time of failure
                    % is large enough to extrude water?
                    thickness = parameters.Ro - (results.Ri-results.z);
                    thickness_at_failure = interp1(results.time,thickness,results.failure_time(1:ifail));
                    overpressure = results.failure_P;% - 1000*parameters.g*thickness_at_failure;
                    
                    for i=1:length(results.failure_eruption_time)
                        % find closest overpressure value
                        [~,ind] = min( abs(results.failure_time - results.failure_eruption_time(i)) );
                        
                        if (rho_i*parameters.g*thickness_at_failure(ind) + overpressure(ind))/(1000*parameters.g) >= thickness_at_failure(ind)
                            plot(results.failure_eruption_time(i)/seconds_in_year/1e6*[1 1],results.failure_erupted_volume(i)/(4*pi*parameters.Ro^2)*[0 1],'b');
                        else
                            plot(results.failure_eruption_time(i)/seconds_in_year/1e6*[1 1],results.failure_erupted_volume(i)/(4*pi*parameters.Ro^2)*[0 1],'b--');
                        end
                        %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_volumechange(i)/(4*pi*Ro^2),'go');
                        %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_pressurechange(i)/(4*pi*Ro^2),'rx');
                    end
                    ylabel('Cracks Reach Ocean');
                    xlabel('Time (Myr)');
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
        end
        % End main code
    end% end loop over moons
end% end postprocess branch
