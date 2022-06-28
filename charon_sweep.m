clear;
close all;
addpath core;

seconds_in_year = 3.1558e7;

thickness = 2e3;
X0 = 0.03;

Xs = linspace(0,0.2,11);
thicknesses = linspace(2e3,20e3,7);

nthick = length(thicknesses);
nX = length(Xs);
results = cell( nthick,nX );
ind=1;
all_p = cell(nthick,nX);
for ithick = 1:nthick
    for iX = 1:nX
        
        p = struct();
        p.Ro = 6.06e5;            % outer radius of ice shell (m)
        p.Ri = p.Ro-thicknesses(ithick);  % (initial) inner radius of ice shell (m)
        p.Rc = p.Ro-2.30e5;         % core radius (m)
        p.g = 0.279;      % used to calculate failure, m/s/s
        p.Ts=40; % Surface temperature (K)
        
        p.Qbelow = @(time) 3e-3; % additional basal heat flux production in W/m^2
        p.relaxation_parameter=1e-2; % used in nonlinear loop.
        p.X0 = Xs(iX);
        p.label = 'Charon';
        p.t_end = 5e3*seconds_in_year;
        all_p{ind} = p;
        ind = ind +1;
    end
end

nrun = nthick*nX;
parfor i=1:nrun 
    results{i} = main_thickening_ice_shell(all_p{i});
end
%%
result = results;
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
xscale = 'log';
figure();
t=tiledlayout(4,1,'TileSpacing','none','Padding','none');
%         t.Units = 'centimeters';
%         t.OuterPosition = [1 1 9.5 9];
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
xlabel('Time (years)');
title(result.label);
ylabel('Depth (km)');
set(gca,'XScale',xscale);
hold on;
for i=1:length(result.failure_time)
    plot(result.failure_time(i)*1e6*[1 1],[result.failure_top(i) result.failure_bottom(i)]/1e3,'r');
end
nexttile
plot(result.time/seconds_in_year,result.Pex);
ylabel('P_{ex} (Pa)');
set(gca,'XScale',xscale);
ax2 = gca();
ax2.Position(3) = ax1.Position(3);
ax2.XLim = ax1.XLim;
ax2.FontSize=8;
hold on
plot(result.failure_time*1e6,result.failure_P,'r.');
end_color = [0 0.9 0];
plot(result.failure_time*1e6,(result.failure_P+result.failure_dP),'LineStyle','none','Color',end_color,'Marker','o','MarkerFaceColor',end_color,'MarkerSize',2);
text(0.025,0.85,char('B'),'FontSize',8,'Units','normalized');
plot(result.failure_time*1e6,result.failure_Pex_crit,'b--');


xlabel('Time (years)');
nexttile
hold on;
for i=1:length(result.failure_time)
    if isnan(result.failure_erupted_volume(i))
        % plot nothing
    else
        if result.failure_P(i) - result.failure_Pex_crit(i) > 0
            plot(result.failure_time(i)*1e6*[1 1],[0 1],'b');
        else
            plot(result.failure_time(i)*1e6*[1 1],[0 1],'b--');
        end
    end
    %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_volumechange(i)/(4*pi*Ro^2),'go');
    %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_pressurechange(i)/(4*pi*Ro^2),'rx');
end
ylabel('Eruption?');
xlabel('Time (years)');
set(gca,'XScale',xscale);
ax3=gca();
ax3.XLim = ax1.XLim;
ax3.Position(3) = ax1.Position(3);
ax3.Box = 'on';
ax3.FontSize=8;
text(0.025,0.85,char('C'),'FontSize',8,'Units','normalized');

nexttile;
plot(result.time/seconds_in_year,result.XNH3,'k');
set(gca,'XScale',xscale);
ylabel('X_{NH_3}')
xlabel('Time (years)');
set(gca,'XLim',ax1.XLim);

fig = gcf();
fig.PaperUnits = 'centimeters';
fig.PaperPosition(3) = 6.00;
fig.Color = 'w';
% filename = sprintf('%s_thickening_nh3-%f_h0-%f.eps',label,initial_ammonia(iAmmonia),...
%     thicknesses(ithick));
% exportgraphics(gcf,filename,'ContentType','vector');

