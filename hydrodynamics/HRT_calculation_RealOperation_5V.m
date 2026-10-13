%% Real-operation Hydraulic Retention Time (HRT) Calculator — 5V
% =========================================================================
% PURPOSE
% Calculates the mean HRT of each vessel in the 5-vessel (5V) gut reactor
% by simulating the ACTUAL pump schedule event by event (feeds, transfers,
% sampling, effluent), instead of approximating it with continuous flow.
%
% OPERATING SCHEDULE MODELLED (confirmed with Fathima, 2026-09-29)
%   - 3 meals/day at 0, 4, 8 h; no pumping from 12 h to 24 h.
%   - Each meal: 60 mL feed + 3 mL pepsin into Stomach at t = 0.
%       +1.00 h  Stomach -> SI  (63 mL); 10 mL pancreatic juice added to SI
%       +1.75 h  SI      -> AC  (73 mL)
%       +2.50 h  AC      -> TC
%       +3.25 h  TC      -> DC
%       +4.00 h  DC      -> effluent (at the same time as the next feed;
%                               for the 3rd meal this is at 12 h)
%   - Sampling: ONCE per day, with the MORNING meal only. 10 mL is
%     withdrawn from AC, TC and DC at the moment each is pumped out.
%       Morning meal : AC->TC 63, TC->DC 53, DC->effluent 43 (+10 sample each)
%       Meals 2 and 3: 73 mL passes straight through AC, TC, DC
%   - Every vessel returns to a 150 mL base volume after each transfer.
%   - Each vessel is perfectly mixed at its current volume when pumped out.
%
% HRT DEFINITIONS REPORTED
%   (1) Per-vessel HRT for a tracer pulse entering with meal 1, 2 or 3:
%       mass-weighted mean exit time (outflow + sampling) of that pulse.
%   (2) Flow-weighted mean HRT: average of (1) weighted by the volume that
%       enters with each meal. This is the true mean residence time of
%       material passing through the vessel.
%   (3) Check: time-averaged volume / time-averaged throughput (Little's
%       law). (2) and (3) must agree; if they don't, something is wrong.
%   (4) Whole-system (cascade) residence time: tracer in the feed, tracked
%       from Stomach to effluent/sampling.
%
% Just run the whole script. Edit only the PARAMETERS section.
% =========================================================================

clear; clc;

%% ================== PARAMETERS — EDIT HERE ==================
P.Vbase      = 150;       % mL, base volume of every vessel
P.feed       = 60;        % mL feed per meal
P.pepsin     = 3;         % mL pepsin added to Stomach with each meal
P.pancreatic = 10;        % mL pancreatic juice added to SI each meal
P.meal_times = [0 4 8];   % h, meal start times within each day

% Pump-out time of each vessel after meal start (h):
%            Stomach  SI    AC    TC    DC
P.t_out   = [1.00    1.75  2.50  3.25  4.00];

% Volume pumped to the NEXT vessel (DC: to effluent), mL.
% Rows = meal 1 (morning), meal 2, meal 3. Columns = Stomach SI AC TC DC
P.out_vol = [63 73 63 53 43;
             63 73 73 73 73;
             63 73 73 73 73];

% Sample volume withdrawn at pump-out, mL (same layout)
P.samp_vol = [0 0 10 10 10;
              0 0  0  0  0;
              0 0  0  0  0];

days  = 60;               % simulation horizon (tracer fully washed out)
names = {'Stomach','SI','AC','TC','DC'};


%% ================== DERIVED QUANTITIES ==================
% Time each vessel receives its inflow, relative to meal start (h)
entry_offset = [0, P.t_out(1:4)];
% How long each vessel stays at its peak volume per meal (h)
dwell = P.t_out - entry_offset;
% Volume entering each vessel per meal (rows = meals, cols = vessels)
in_vol = [(P.feed + P.pepsin)*ones(3,1), ...
          P.out_vol(:,1) + P.pancreatic, ...
          P.out_vol(:,2:4)];

% Mass-balance check: what enters each vessel per meal must leave
bal = in_vol - (P.out_vol + P.samp_vol);
if any(abs(bal(:)) > 1e-9)
    error('Volume balance fails: inflow ~= outflow + sample for some vessel/meal.');
end


%% ================== (1)+(2) PER-VESSEL HRT ==================
HRT_meal  = zeros(3,5);
recovered = zeros(3,5);
for i = 1:5
    for k = 1:3
        t0 = P.meal_times(k) + entry_offset(i);
        [te, me] = simulate_tracer(P, days, t0, i, true);
        recovered(k,i) = sum(me);
        HRT_meal(k,i)  = sum(me .* (te - t0)) / sum(me);
    end
end
HRT_fw = sum(HRT_meal .* in_vol, 1) ./ sum(in_vol, 1);


%% ================== (3) LITTLE'S LAW CHECK ==================
Vmean  = P.Vbase + sum(in_vol, 1) .* dwell / 24;   % time-averaged volume (mL)
Q      = sum(in_vol, 1) / 24;                      % mean throughput (mL/h)
HRT_LL = Vmean ./ Q;


%% ================== (4) WHOLE-SYSTEM CASCADE ==================
HRT_sys = zeros(1,3);
for k = 1:3
    t0 = P.meal_times(k);
    [te, me] = simulate_tracer(P, days, t0, 1, false);
    HRT_sys(k) = sum(me .* (te - t0)) / sum(me);
end


%% ================== OUTPUT ==================
fprintf('=== 5V REAL-OPERATION HRT (h) ===\n');
fprintf('%-8s %8s %8s %8s | %10s %10s | %8s\n', 'Vessel', 'Meal 1', 'Meal 2', 'Meal 3', ...
        'Flow-wtd', 'Little V/Q', 'Vmean');
for i = 1:5
    fprintf('%-8s %8.2f %8.2f %8.2f | %10.2f %10.2f | %8.2f\n', names{i}, ...
        HRT_meal(1,i), HRT_meal(2,i), HRT_meal(3,i), HRT_fw(i), HRT_LL(i), Vmean(i));
end
fprintf('%-8s %8.2f %8.2f %8.2f | %10.2f %10.2f\n', 'SUM', sum(HRT_meal,2), ...
        sum(HRT_fw), sum(HRT_LL));
fprintf('\nWhole system (feed -> effluent/sampling), tracer in meal 1/2/3: %.2f / %.2f / %.2f h\n', HRT_sys);
fprintf('Mean over the three meals: %.2f h\n', mean(HRT_sys));
fprintf('\nTracer recovery (should all be 1.0000): min = %.6f, max = %.6f\n', ...
        min(recovered(:)), max(recovered(:)));
fprintf('Max |flow-weighted - Little| = %.2e h (should be ~0)\n', max(abs(HRT_fw - HRT_LL)));

% Previously reported (Chapter 3 / Appendix 1), for comparison:
%   Stomach 16.9 h | AC 14.6 h | TC 16.9 h | DC 20.1 h


%% ================== ENGINE ==================
function ev = build_events(P, days)
    % Rows: [time, type (1 = feed, 2 = pump-out), meal index, vessel]
    ev = zeros(days*3*6, 4);
    r = 0;
    for d = 0:days-1
        for k = 1:3
            T = 24*d + P.meal_times(k);
            r = r + 1; ev(r,:) = [T, 1, k, 1];
            for i = 1:5
                r = r + 1; ev(r,:) = [T + P.t_out(i), 2, k, i];
            end
        end
    end
    ev = sortrows(ev, [1 2]);
end

function [t_exit, m_exit] = simulate_tracer(P, days, inj_time, inj_vessel, isolate)
    % Unit tracer is added to vessel inj_vessel at the moment it receives
    % its inflow at time inj_time. If isolate = true, tracer leaving that
    % vessel counts as exited (per-vessel HRT). If false, it passes on down
    % the chain (whole-system HRT). Sampled tracer always counts as exited.
    ev = build_events(P, days);
    V = P.Vbase * ones(1,5);
    M = zeros(1,5);
    t_exit = []; m_exit = [];
    injected = false;

    for e = 1:size(ev,1)
        t = ev(e,1); type = ev(e,2); k = ev(e,3); i = ev(e,4);

        if type == 1
            % Feed + pepsin into Stomach
            V(1) = V(1) + P.feed + P.pepsin;
            if ~injected && inj_vessel == 1 && abs(t - inj_time) < 1e-9
                M(1) = M(1) + 1; injected = true;
            end
        else
            % Pump-out of vessel i (perfectly mixed at current volume)
            out = P.out_vol(k,i); smp = P.samp_vol(k,i); tot = out + smp;
            m = M(i) * tot / V(i);
            M(i) = M(i) - m;  V(i) = V(i) - tot;
            m_out = m * out / tot;  m_smp = m - m_out;

            if smp > 0
                t_exit(end+1) = t; m_exit(end+1) = m_smp; %#ok<AGROW>
            end
            if i == 5 || (isolate && i == inj_vessel)
                t_exit(end+1) = t; m_exit(end+1) = m_out; %#ok<AGROW>
            else
                M(i+1) = M(i+1) + m_out;
            end

            if i < 5
                V(i+1) = V(i+1) + out + (i == 1) * P.pancreatic;
                if ~injected && inj_vessel == i+1 && abs(t - inj_time) < 1e-9
                    M(i+1) = M(i+1) + 1; injected = true;
                end
            end
        end
    end
end
