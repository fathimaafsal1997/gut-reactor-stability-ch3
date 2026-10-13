%% Washout curves per vessel: nominal (continuous CSTR) vs effective (real pulse feeding) — 5V
% =========================================================================
% For each vessel on its own, plots the % of food still in the vessel
% against time since it entered:
%   - dashed grey : ideal continuous CSTR at the NOMINAL HRT
%                   remaining = exp(-t / HRT_nom)
%   - solid blue  : real pump schedule (pulse feeding, overnight pause),
%                   averaged over the three meals of a day, weighted by
%                   the volume each meal brings in
% The area under each curve equals that curve's HRT:
%   area under the dashed curve = HRT_nom
%   area under the blue curve   = HRT_eff (= HRT_nom + tau)
% The script prints both, so you can check this.
% =========================================================================
clear; clc; close all;

%% ================== PARAMETERS ==================
Vb    = 150;                              % base volume (mL)
meals = [0 4 8];                          % meal times (h)
%          Stomach  SI    AC    TC    DC
t_in  = [0.00     1.00  1.75  2.50  3.25];  % when each vessel receives its inflow (h after meal start)
t_out = [1.00     1.75  2.50  3.25  4.00];  % when each vessel pumps out
% Volume entering each vessel per meal (rows = meal 1,2,3).
% Morning sampling (10 mL from AC and TC) reduces what reaches TC and DC at meal 1.
inflow = [63 73 73 63 53;
          63 73 73 73 73;
          63 73 73 73 73];
names = {'Stomach','Small intestine','Ascending colon','Transverse colon','Descending colon'};

tgrid = (0:0.02:120)';                    % plotting grid (h since entry)
days  = 200;                              % long horizon for the area check

%% ================== COMPUTE ==================
fig = figure('Color','w','Units','centimeters','Position',[1 1 38 23]);
try, theme(fig, 'light'); end   % force light theme (MATLAB R2025a+ dark mode)
tl  = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
panel = 'ABCDE';

for i = 1:5
    w   = inflow(:,i);                    % weights = volume entering with each meal
    Q   = sum(w) / 24;                    % mean flow (mL/h)
    nom = Vb / Q;                         % nominal HRT
    tau = t_out(i) - t_in(i);
    eff = nom + tau;                      % effective HRT

    Rgrid = zeros(size(tgrid)); area_real = 0;
    for k = 1:3
        [ts, rem] = staircase(i, k, Vb, meals, t_in, t_out, inflow, days);
        idx   = discretize(tgrid, [ts; inf]);        % which step each grid time falls in
        Rgrid = Rgrid + w(k) * rem(idx);
        area_real = area_real + w(k) * sum(rem(1:end-1) .* diff(ts));
    end
    Rgrid = Rgrid / sum(w);
    area_real = area_real / sum(w);
    fprintf('%-17s nominal %.2f h | effective %.2f h | area under real curve %.2f h\n', ...
        names{i}, nom, eff, area_real);

    ax = nexttile; hold on
    plot(tgrid, 100*exp(-tgrid/nom), '--', 'Color',[0.48 0.48 0.48], 'LineWidth',2, ...
        'DisplayName','Continuous CSTR (nominal)');
    plot(tgrid, 100*Rgrid, '-', 'Color',[0.165 0.471 0.839], 'LineWidth',2, ...
        'DisplayName','Real pulse feeding (effective)');
    xline(nom, '--', 'Color',[0.48 0.48 0.48], 'LineWidth',1.2, 'HandleVisibility','off');   % nominal HRT
    xline(eff, '-',  'Color',[0.165 0.471 0.839], 'LineWidth',1.2, 'HandleVisibility','off'); % effective HRT
    text(0.97, 0.96, sprintf('Nominal %.2f h\nEffective %.2f h', nom, eff), 'Units','normalized', ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontWeight','bold', 'FontSize',12);
    title([panel(i) '  ' names{i}], 'FontSize',14, 'HorizontalAlignment','left', ...
        'Units','normalized', 'Position',[0 1.01]);
    xlim([0 96]); ylim([0 102]); xticks(0:24:96);
    if i >= 3, xlabel('Time since entry (h)'); end
    if i == 1 || i == 4, ylabel('Food remaining in vessel (%)'); end
    ax.FontSize = 13; ax.FontWeight = 'bold'; ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 1.2; ax.Box = 'off'; ax.YGrid = 'on'; ax.GridColor = [0.85 0.85 0.85];
end

% ---------- 6th tile: legend + note ----------
ax6 = nexttile; axis(ax6, 'off'); hold(ax6, 'on')
p1 = plot(ax6, NaN, NaN, '--', 'Color',[0.48 0.48 0.48], 'LineWidth',2);
p2 = plot(ax6, NaN, NaN, '-',  'Color',[0.165 0.471 0.839], 'LineWidth',2);
legend(ax6, [p1 p2], {'Continuous CSTR (nominal HRT)', 'Real pulse feeding (effective HRT)'}, ...
    'Location','north', 'Box','off', 'FontWeight','bold', 'FontSize',14);
text(ax6, 0.5, 0.25, {'Vertical lines mark each HRT;', 'area under each curve = its HRT'}, ...
    'Units','normalized', 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',12);

exportgraphics(fig, 'HRT_washout_nominal_vs_effective_5V.png', 'Resolution', 300, 'BackgroundColor', 'white');

%% ================== FUNCTION ==================
function [ts, rem] = staircase(i, k, Vb, meals, t_in, t_out, inflow, days)
    % Food entering vessel i with meal k (vessel on its own).
    % At each pump-out the vessel holds Vb + inflow and pumps out 'inflow'
    % (to the next vessel + sample), so a fraction Vb/(Vb + inflow) stays.
    % Returns step times (h since entry) and the fraction remaining after each step.
    t0 = meals(k) + t_in(i);
    T  = []; K = [];
    for d = 0:days-1
        for kk = 1:3
            T(end+1) = 24*d + meals(kk) + t_out(i); %#ok<AGROW>
            K(end+1) = kk;                          %#ok<AGROW>
        end
    end
    keep = T > t0 + 1e-9;  T = T(keep); K = K(keep);
    [T, order] = sort(T);  K = K(order);
    ts = [0; (T(:) - t0)];
    r  = Vb ./ (Vb + inflow(K, i));                % fraction staying at each pump-out
    rem = [1; cumprod(r(:))];
end
