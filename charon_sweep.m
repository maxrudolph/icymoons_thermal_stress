clear;
close all;
addpath core;

seconds_in_year = 3.1558e7;

% Xs = linspace(0.03,0.1,2);% 0-0.2, 7 works
Xs = [0 0.03 0.1]

thicknesses = linspace(2e3,2e3,1);% 5 works
% loop for initial thickness
nthick = length(thicknesses);
nX = length(Xs);
% loop for surface temperature
nTs = 2;
Ts = [40 60];
% loop for melting point viscosity
mubs = [1e13 1e14 1e15];
nmub = 3;

results = cell( nthick*nX*nTs*nmub,1 );
ind=1;
all_p = cell(nthick*nX*nTs*nmub,1);
irun = 1;
for iTs = 1:nTs
    for ithick = 1:nthick
        for iX = 1:nX      
            for imub = 1:nmub
                p = struct();
                p.Ro = 6.06e5;            % outer radius of ice shell (m)
                p.Ri = p.Ro-thicknesses(ithick);  % (initial) inner radius of ice shell (m)
                p.Rc = p.Ro-2.30e5;         % core radius (m)
                p.g = 0.279;      % used to calculate failure, m/s/s
                p.Ts=Ts(iTs); % Surface temperature (K)
                p.mub = mubs(imub);
                p.Qbelow = @(time) 0*3e-3; % additional basal heat flux production in W/m^2
                p.relaxation_parameter=1e-2; % used in nonlinear loop.
                p.X0 = Xs(iX);
                p.label = 'Charon';
                p.t_end = 5e8*seconds_in_year;
                all_p{ithick,iX} = p;
                irun = irun+1;
            end
        end
    end
end

nrun = nthick*nX;
parfor i=1:nrun
    results{i} = main_thickening_ice_shell(all_p{i});
end

save('charon_08042022.mat','-v7.3');

%% summary plots
all_cracks_reach_ocean = zeros(size(results));
all_max_failure_thickness = zeros(size(results));
all_X0 = zeros(size(results));

nthick = size(results,1);
nX = size(results,2);

for i=1:length(results(:))
    result = results{i};
    all_X0(i) = result.parameters.X0;
    if isempty(result.failure_erupted_volume) || all(isnan(result.failure_erupted_volume))
        all_cracks_reach_ocean(i) = 0;
        all_max_failure_thickness(i) = NaN;
    else
        all_cracks_reach_ocean(i) = 1;
        all_max_failure_thickness(i) = max(result.failure_z( ~isnan(result.failure_erupted_volume) ) );
    end
end
figure();
imagesc(Xs,thicknesses,all_cracks_reach_ocean);
colorbar();
caxis([0 1]);
xlabel('Initial ammonia content (X_0)');
ylabel('Initial thickness (km)');
set(gca,'YDir','normal');

figure();
imagesc(Xs,thicknesses,all_max_failure_thickness/1e3);
hcb=colorbar();
hcb.Label.String = 'Thickest ice when cracks reach ocean (km)';
caxis([0 15])
xlabel('Initial ammonia content (X_0)');
ylabel('Initial thickness (km)');
set(gca,'YDir','normal');
set(gcf,'Color','w');
exportgraphics(gcf,'Charon-thickest-cracks.eps','ContentType','Vector');


figure();
imagesc(Xs,thicknesses,all_max_failure_thickness);%,10);
colorbar();
caxis([0 15000])
xlabel('Initial ammonia content (X_0)');
ylabel('Initial thickness (km)');
set(gca,'YDir','normal');


%%
result = results{1,1};
figure();

subplot(3,1,1);
plot(result.time/seconds_in_year,result.z)
ylabel('Amount of freezing (m)');
subplot(3,1,2);
plot(result.time/seconds_in_year,result.sigma_t(1,:),'DisplayName',sprintf('%.02f km',result.save_depths(1)/1000));
hold on
plot(result.time/seconds_in_year,result.sigma_t(10,:),'DisplayName',sprintf('%.02f km',result.save_depths(10)/1000));
plot(result.time/seconds_in_year,result.sigma_t(20,:),'DisplayName',sprintf('%.02f km',result.save_depths(20)/1000));
legend();
ylabel('\sigma_t (Pa)')
subplot(3,1,3);
plot(result.time/seconds_in_year,result.Pex);
ylabel('Overpressure (Pa)');

%% Pseudocolor stress plot
for ithick = [1 nthick]
    for iX = [1 nX]
        
        result = results{ithick,iX}
        xscale = 'log';
        xtick = [1e3 1e4 1e5 1e6 1e7 1e8];
%         xscale = 'linear'
%         xtick = 'auto'

        figure();
        f=gcf();
        f.Position(3:4) = [945 890];
        t=tiledlayout(4,1,'TileSpacing','tight','Padding','none');
        t.Units = 'centimeters';
        t.OuterPosition = [1 1 10 14];
        nexttile
        
        contourf(result.time/seconds_in_year,result.save_depths/1000,result.sigma_t,64,'Color','none'); %shading flat;
        hold on
        plot(result.time/seconds_in_year,((result.parameters.Ro-result.parameters.Ri)+result.z)/1000,'Color','k','LineWidth',1);
        %         set(gca,'YLim',[0 ceil(1+max(((Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        ax1 = gca();
        ax1.FontSize=8;
        hcb = colorbar();
        hcb.Label.String = 'Tensile Stress (Pa)';
        text(0.025,0.85,char('A'),'FontSize',8,'Units','normalized');
        %         xlabel('Time (years)');
        title(result.label);
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        %         set(gca,'XTick',xtick);
        hold on;
        for i=1:length(result.failure_time)
            plot(result.failure_time(i)*1e6*[1 1],[result.failure_top(i) result.failure_bottom(i)]/1e3,'r');
        end
        nexttile
        plot(result.time/seconds_in_year,result.Pex/1e6);
        ylabel('P_{ex} (MPa)');
        set(gca,'XScale',xscale);
        ax2 = gca();
        %         ax2.Position(3) = ax1.Position(3);
        ax2.XLim = ax1.XLim;
        ax2.FontSize=8;
        hold on
        plot(result.failure_time*1e6,result.failure_P/1e6,'r.');
        end_color = [0 0.9 0];
        plot(result.failure_time*1e6,(result.failure_P/1e6+result.failure_dP/1e6),'LineStyle','none','Color',end_color,'Marker','o','MarkerFaceColor',end_color,'MarkerSize',2);
        text(0.025,0.85,char('B'),'FontSize',8,'Units','normalized');
        plot(result.failure_time*1e6,result.failure_Pex_crit/1e6,'k');
        
        
        %         xlabel('Time (years)');
%         nexttile
        hold on;
        for i=1:length(result.failure_time)
            if isnan(result.failure_erupted_volume(i))
                % plot nothing
            else
                if result.failure_P(i) - result.failure_Pex_crit(i) > 0
                    plot(result.failure_time(i)*1e6*[1 1],ax2.YLim,'b');
                else
                    plot(result.failure_time(i)*1e6*[1 1],ax2.YLim,'b--');
                end
            end
            %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_volumechange(i)/(4*pi*Ro^2),'go');
            %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_pressurechange(i)/(4*pi*Ro^2),'rx');
        end
%         ylabel('Eruption?');
        %         xlabel('Time (years)');
        set(gca,'XScale',xscale);
        ax3=gca();
        ax3.XLim = ax1.XLim;
        %         ax3.Position(3) = ax1.Position(3);
        ax3.Box = 'on';
        ax3.FontSize=8;
%         text(0.025,0.85,char('C'),'FontSize',8,'Units','normalized');
        
        nexttile;
        plot(result.time/seconds_in_year,result.XNH3(1:length(result.time)),'k');
        set(gca,'XScale',xscale);
        ylabel('X (NH_3)')
        xlabel('Time (years)','FontSize',8);
        
        set(gca,'XLim',ax1.XLim);
        set(gca,'XTick',ax1.XTick);
        
        text(0.025,0.85,char('C'),'FontSize',8,'Units','normalized');
        
        ax4=gca();
        ax4.FontSize=8;
        fig = gcf();
        %fig.PaperUnits = 'centimeters';
        %fig.PaperPosition(3) = 6.00;
        fig.Color = 'w';
        filename = sprintf('%s_thickening_nh3-%f_h0-%f.eps',result.label,result.parameters.X0,...
            result.parameters.Ro-result.parameters.Ri);
        exportgraphics(t,filename,'ContentType','vector');
    end
end
