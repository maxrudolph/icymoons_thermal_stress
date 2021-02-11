%% This function sweeps parameter space, running a series of models
%Europa-like conditions
clear;
close all;
addpath core;
seconds_in_year = 3.1558e7;
do_runs = false
if do_runs
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
            parameters.k = @(T) 651./T; % Petrenko and Whitworth 1999
            parameters.label = 'Europa';
        else
            % Enceladus
            parameters = struct();
            parameters.Tb = 273;
            parameters.Ts = 100;
            parameters.k = @(T) 651./T; %Petrenko and Whitworth 1999

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
        % parfor here for real runs:
        parfor irun=1:ndQ*nthick
            all_results{irun} = main_cyclic_thermomechanical_model(all_parameters{irun});
        end
        save([parameters.label '_workspace.mat'],'all_parameters','all_results','ndQ','nthick','thicknesses','dQ','-v7.3');
    end
else% postprocess:
    %% Load results
    close all;
    addpath ~/sw/matlab/crameri

    for moon=0:1
        clearvars -except seconds_in_year moon;
        seconds_in_year = 3.1558e7;
        if moon==0
            load('~/Europa_workspace.mat');
        else
            load('~/Enceladus_workspace.mat');
        end
        
        
        %% Analyze models
        all_erupted_volume = zeros(size(all_results));
        all_failure_events = zeros(size(all_results));
        all_failure_fraction = zeros(size(all_results));
        all_actual_thickness = zeros(size(all_results));
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
                all_erupted_volume(idQ,ithick) = sum( results.failure_erupted_volume(time_mask) );
                all_failure_events(idQ,ithick) = nnz( results.failure_time(time_mask) );
                % calculate the average thickness.
                thickness = parameters.Ro-(results.Ri-results.z);
                average_thickness = 1/(results.time(end)-results.time(1))*sum( (thickness(1:end-1) + thickness(2:end))/2 .* diff(results.time));
                all_actual_thickness(idQ,ithick) = average_thickness;
                % approximate ice shell thickness at time of failure
                if any(results.failure_thickness)
                    failure_thickness = interp1(results.time,parameters.Ro-(results.Ri-results.z),results.failure_time(time_mask));
                    all_failure_fraction(idQ,ithick) = max(results.failure_thickness ./ failure_thickness);
                end
            end
        end
        if moon==0
            figure(901);
            clf;
            t1=tiledlayout(3,2,'Padding','none','TileSpacing','none','Units','centimeters');
            t1.OuterPosition(3:4) = [18 12];
        else
            figure(901);
        end
        
   
        nexttile(1+moon);
        dt = max(results.time) - 2*parameters.perturbation_period;
        contourf(thicknesses/1e3,dQ,all_erupted_volume/(4*pi*parameters.Ro^2)/(dt/seconds_in_year/1e6),16,'Color','none');
        hcb=colorbar();
        hcb.Label.String = 'Erupted Volume (m/Myr)';
        colormap(crameri('davos'));
        ylabel('\Delta q/q_0');
        xlabel('Equilibrium Thickness (km)');
        if( max(all_erupted_volume(:))==0 )
            caxis([0 1])
        end
        title(parameters.label);
        
        % Panel 2
        nexttile(3+moon);
        contourf(thicknesses/1e3,dQ,all_failure_fraction,16,'Color','none');
        caxis([0 1]);
        hcb = colorbar();
        hcb.Label.String = 'Fractional penetration';
        ylabel('\Delta q/q_0');
        xlabel('Equilibrium Thickness (km)');
        
        % Panel 3
        nexttile(5+moon);
        dt = max(results.time) - 2*parameters.perturbation_period;
        contourf(thicknesses/1e3,dQ,all_failure_events/(dt/parameters.perturbation_period),0:20,'Color','none');
        hcb=colorbar();
        
        hcb.Label.String = 'Failure events per cycle';
        ylabel('\Delta q/q_0');
        xlabel('Equilibrium Thickness (km)');
        colormap(crameri('davos'));

        if moon == 1
            exportgraphics(gcf,[parameters.label '_regimes.eps']);
        end
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
%                 ss_thickness = k*(273-100)./results.qb;
               
                qsurf = results.Qtot./(4*pi*parameters.Ro^2);
                ss_thickness = find_steady_h(parameters.Ro,parameters.Tb,parameters.Ts,qsurf);
                if any( imag(ss_thickness) ~= 0 )
                   ss_thickness = real(ss_thickness);
                   warning('complex thickness encountered - this case is not reasonable');
                end
                [ss_peaks,ss_ind] = findpeaks(ss_thickness,'MinPeakProminence',std(ss_thickness)/2);
                [actual_peaks,actual_ind] = findpeaks(actual_thickness,'MinPeakProminence',std(actual_thickness)/2);
                phase_lag = median(results.time(actual_ind) - results.time(ss_ind))/seconds_in_year/1e6;
                all_phase_lags(idQ,ithick) = phase_lag;
                all_ampfrac(idQ,ithick) = (max(actual_thickness)-min(actual_thickness))/(max(ss_thickness)-min(ss_thickness));
                %         plot(results.time/seconds_in_year,ss_thickness,'--');
            end
        end
        figure;
        t=tiledlayout(2,1,'Padding','none','TileSpacing','none','Units','centimeters');
        t.OuterPosition(3:4) = [9 9];
        
        nexttile
        contourf(thicknesses/1e3,dQ,all_phase_lags,20,'Color','none'); crameri('davos')
        hcb=colorbar();
        hcb.Label.String='Lag (Myr)';
        ylabel('\Delta q/q_0');
        xlabel('Average Thickness (km)');
        title(parameters.label);
        
        nexttile
        contourf(thicknesses/1e3,dQ,all_ampfrac,20,'Color','none'); crameri('davos')
        hcb=colorbar();
        hcb.Label.String='\Delta h/\Delta h_{eq}';
        ylabel('\Delta q/q_0');
        xlabel('Average Thickness (km)');
        caxis([0 1])
        exportgraphics(gcf,[parameters.label '_lag_eqlm.eps']);
        
        %% Plot maximum fractional penetration
        figure();
        
        
        %% Plot outcomes of individual runs
        for ithick=1:1
            for idQ=1:3:ndQ
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
    end% end loop over moons
end% end postprocess branch