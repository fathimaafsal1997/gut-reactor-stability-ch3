%% Feeding and starvation cycle with the HRT variables marked (5V, real operation)
% =========================================================================
% Same step plot as HRT_feeding_starvation_schedule.m, with the variables of
% the HRT derivation drawn on it:
%   V        base volume of each compartment (150 mL)
%   V_flow   volume transferred at each pump-out (step height)
%   T_peak   time the inflow sits in the vessel before it is pumped on (pulse width)
%   t1,t2,t3 gaps between successive transfers (4, 4, 16 h; same for every vessel)
%   T        length of one feeding cycle (24 h)
%   Q_av     mean flow = sum of V_flow over one cycle / T
% Stomach panel (day 1) carries the arrows; every panel lists its values.
% Uses only stairs / patch / line / text so it draws on all MATLAB versions.
% =========================================================================
clear; clc; close all;

%% ================== PARAMETERS ==================
Vb    = 150;                       % V, base volume (mL)
meals = [0 4 8];                   % feed times within the day (h)
T     = 24;                        % cycle length (h)
days  = 2;  tend = 48;
%          Stomach  SI    AC    TC    DC
t_in  = [0.00     1.00  1.75  2.50  3.25];   % receives inflow (h after feed)
t_out = [1.00     1.75  2.50  3.25  4.00];   % transfers out
inflow = [63 73 73 63 53;          % V_flow, feed 1 (morning; 10 mL samples taken from AC, TC, DC)
          63 73 73 73 73;          % feed 2
          63 73 73 73 73];         % feed 3
names = {'Stomach','Small intestine','Ascending colon','Transverse colon','Descending colon'};
cols  = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.929 0.631 0.000; 0.910 0.482 0.643];
k     = [0.85 0.15 0.15];          % annotation colour (dark red)

Qav   = sum(inflow,1) / T;         % mL/h
Tpeak = t_out - t_in;              % h

%% ================== FIGURE ==================
fig = figure('Color','w','Units','centimeters','Position',[2 1 34 31]);
try, theme(fig, 'light'); end      % force light theme (MATLAB R2025a+ dark mode)
tl = tiledlayout(6,1,'TileSpacing','compact','Padding','compact');

for i = 1:5
    if i == 1, ax = nexttile([2 1]); yl = [140 300]; else, ax = nexttile; yl = [140 235]; end
    hold(ax,'on');
    % starvation shading
    for s = [12 36]
        patch(ax, [s s+12 s+12 s], [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], 'EdgeColor','none');
    end
    % volume step series
    ev = [];
    for d = 0:days-1
        for m = 1:3
            Tm = 24*d + meals(m);
            ev = [ev; Tm + t_in(i), +inflow(m,i); Tm + t_out(i), -inflow(m,i)]; %#ok<AGROW>
        end
    end
    ev = sortrows(ev, 1);
    ts = [0; ev(:,1); tend];
    vs = Vb + [0; cumsum(ev(:,2)); sum(ev(:,2))];
    stairs(ax, ts, vs, '-', 'Color', cols(i,:), 'LineWidth', 2.2);

    % name (day-1 band) and values (day-2 band)
    text(ax, 18, 192, names{i}, 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',13);
    if all(inflow(:,i) == inflow(1,i))
        vf = sprintf('%d', inflow(1,i));
    else
        vf = sprintf('%d/%d/%d', inflow(:,i));
    end
    text(ax, 42, 200, sprintf('V = %d mL,  V_{flow} = %s mL', Vb, vf), ...
        'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);
    text(ax, 42, 178, sprintf('Q_{av} = %.2f mL/h,  T_{peak} = %.2f h', Qav(i), Tpeak(i)), ...
        'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);

    % sample markers (morning feed only) for AC, TC, DC
    if i >= 3
        for d = 0:days-1
            text(ax, 24*d + t_out(i), 228, char(9660), 'HorizontalAlignment','center', 'FontSize',14);
        end
    end

    % ---------- stomach panel: variables drawn on day 1 ----------
    if i == 1
        for d = 0:days-1
            for m = 1:3
                text(ax, 24*d + meals(m), 220, 'feed \downarrow', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', 'FontWeight','bold', 'FontSize',10);
            end
        end
        text(ax, 42, 225, 'starvation (no pumping)', 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',11);
        text(ax, 18, 225, 'starvation (no pumping)', 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',11);

        % t1, t2, t3 between transfers (pump-outs at 1, 5, 9, 25 h)
        P = [1 5 9 25]; yt = 255;
        for x = P, line(ax, [x x], [230 yt], 'Color',[0.55 0.55 0.55], 'LineStyle',':', 'LineWidth',1.2); end
        lab = {'t_1 = 4 h','t_2 = 4 h','t_3 = 16 h'};
        for j = 1:3
            harrow(ax, P(j), P(j+1), yt, k);
            text(ax, mean(P(j:j+1)), yt + 3, lab{j}, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);
        end
        % T, one cycle (feed to next morning feed)
        yT = 282;
        line(ax, [24 24], [236 yT], 'Color',[0.55 0.55 0.55], 'LineStyle',':', 'LineWidth',1.2);
        harrow(ax, 0, 24, yT, k);
        text(ax, 12, yT + 3, 'T = 24 h (one feeding cycle)', 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);
        % V_flow: step height (second feed)
        varrow(ax, 4.5, Vb, Vb + inflow(2,1), k);
        text(ax, 5.2, 182, 'V_{flow}', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);
        % T_peak: pulse width (third feed)
        harrow(ax, 8, 9, 190, k);
        text(ax, 9.3, 190, 'T_{peak}', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);
        % V: base volume
        text(ax, 1.4, 153, 'V (base)', 'VerticalAlignment','bottom', 'FontWeight','bold', 'FontSize',11.5, 'Color',k);
    end

    xlim(ax, [0 tend]); ylim(ax, yl); yticks(ax, [150 200]); xticks(ax, 0:4:tend);
    if i < 5, ax.XTickLabel = []; end
    ax.FontSize = 13; ax.FontWeight = 'bold'; ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 1.2; ax.Box = 'off'; ax.Layer = 'top'; ax.TickDir = 'out';
    ax.YGrid = 'on'; ax.GridColor = [0.85 0.85 0.85];
end
xlabel(tl, 'Time (h)', 'FontWeight','bold', 'FontSize',13);
ylabel(tl, 'Volume in compartment (mL)', 'FontWeight','bold', 'FontSize',13);
title(tl, 'Feeding and starvation cycle: compartment volumes under real operation', ...
    'FontWeight','bold', 'FontSize',15);
subtitle(tl, {[char(9660) ' 10 mL sample (morning feed only)      grey: starvation, no pumping (12–24 h)'], ...
    'Q_{av} = \SigmaV_{flow} / T      HRT_{nom} = V / Q_{av}      HRT_{eff} = HRT_{nom} + T_{peak}      (t_1, t_2, t_3 and T are the same for every compartment)'}, ...
    'FontWeight','bold', 'FontSize',12);

exportgraphics(fig, 'HRT_feeding_starvation_schedule_annotated.png', 'Resolution', 300, 'BackgroundColor', 'white');

%% ================== FUNCTIONS ==================
function harrow(ax, x1, x2, y, c)          % horizontal double-headed arrow
    line(ax, [x1 x2], [y y], 'Color',c, 'LineWidth',1.6);
    line(ax, x1, y, 'Marker','<', 'MarkerSize',6, 'MarkerFaceColor',c, 'MarkerEdgeColor',c, 'LineStyle','none');
    line(ax, x2, y, 'Marker','>', 'MarkerSize',6, 'MarkerFaceColor',c, 'MarkerEdgeColor',c, 'LineStyle','none');
end

function varrow(ax, x, y1, y2, c)          % vertical double-headed arrow
    line(ax, [x x], [y1 y2], 'Color',c, 'LineWidth',1.6);
    line(ax, x, y1, 'Marker','v', 'MarkerSize',6, 'MarkerFaceColor',c, 'MarkerEdgeColor',c, 'LineStyle','none');
    line(ax, x, y2, 'Marker','^', 'MarkerSize',6, 'MarkerFaceColor',c, 'MarkerEdgeColor',c, 'LineStyle','none');
end
