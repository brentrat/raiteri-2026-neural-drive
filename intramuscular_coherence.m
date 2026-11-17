function [auc_delta, auc_alpha, auc_beta, auc_gamma] = intramuscular_coherence(all_firing_times, cst_subset, fsamp, first)
%% Intramuscular coherence analysis of spike train data across multiple trials
%
%% Details:
%   all_firing_times: 3D cell matrix [num_trials x condition of interest (pass 1 only) x num_units]
%                     Each cell contains a 1D vector of spike sample indices.
%   cst_subset:       CST matrix [num_trials x num_samples] specific to this condition.
%   fsamp:            Sampling frequency (Hz).
%   first:            First cropping index of spike sample indices.

%% Inputs
num_permutations = 1000;          % 100 is common, but 1000 leads to more stable results
win_len = 1 * fsamp;              % 1.0-second Hanning window for 1 Hz resolution
noverlap = round(0.75 * win_len); % 75% overlap
delta_band = [0.5, 5];            % Delta band range (Hz)

%% Data bounds
[num_trials, num_samples] = size(cst_subset);
total_mu = size(all_firing_times, 3);

min_mu_per_pool = 2;
max_mu_per_pool = floor(total_mu / 2);

if max_mu_per_pool < min_mu_per_pool
    error('Dataset only has %d matched motor units. Need at least 4 to create pools.', total_mu);
end

% Determine number of loops through sub-pool sizes
num_sizes = max_mu_per_pool - min_mu_per_pool + 1;
n_freqs = (win_len / 2) + 1;
h_win = hann(win_len);

%% Preallocation
% Pre-calculate the frequency vector 'f'
[~, f] = cpsd(zeros(num_samples, 1), zeros(num_samples, 1), h_win, noverlap, win_len, fsamp);

% Preallocate output variables across sub-pool loops if looping through
% sub-pools (currently not doing this so this code is commented out)
%auc_delta  = zeros(num_sizes, 1);
%auc_alpha  = zeros(num_sizes, 1);
%auc_beta   = zeros(num_sizes, 1);
%auc_gamma  = zeros(num_sizes, 1);

%color_idx = 1;

% Loop through sub-pools
for n_mu = max_mu_per_pool

    % Preallocate matrix for permutations
    perm_coherence = zeros(n_freqs, num_permutations);

    % Run permutations
    for p = 1:num_permutations

        % Preallocate spectral accumulators across trials
        Pxx_total = zeros(n_freqs, 1);
        Pyy_total = zeros(n_freqs, 1);
        Pxy_total = zeros(n_freqs, 1);

        % Loop through trials
        for t = 1:num_trials

            % Find which units actually fired in THIS specific trial and condition
            active_units_in_trial = [];
            for mu = 1:total_mu
                if ~isempty(all_firing_times{t, mu})
                    active_units_in_trial = [active_units_in_trial, mu];
                end
            end

            % Shuffle and select pools strictly from units active in this trial
            shuffled_idx = active_units_in_trial(randperm(length(active_units_in_trial)));

            % Dynamically limit the pool size to what this specific trial can provide
            current_n_mu = min(n_mu, floor(length(active_units_in_trial) / 2));
            
            pool1_idx = shuffled_idx(1:current_n_mu);
            pool2_idx = shuffled_idx(current_n_mu + 1:2*current_n_mu);

            % Allocate clean arrays for the current trial
            cst1_raw = zeros(num_samples, 1);
            cst2_raw = zeros(num_samples, 1);

            % Reconstruct pool 1
            for idx = pool1_idx
                spikes = all_firing_times{t, idx};
                spikes = spikes-first(t, 1)+1;
                cst1_raw(spikes) = cst1_raw(spikes) + 1;
            end

            % Reconstruct pool 2
            for idx = pool2_idx
                spikes = all_firing_times{t, idx};
                spikes = spikes-first(t, 1)+1;
                cst2_raw(spikes) = cst2_raw(spikes) + 1;
            end

            % Detrend individual raw trials before spectral transforms
            cst1_detrend = detrend(cst1_raw);
            cst2_detrend = detrend(cst2_raw);

            % Calculate cross- and auto-spectral densities per trial
            [Pxy_trial, ~] = cpsd(cst1_detrend, cst2_detrend, h_win, noverlap, win_len, fsamp);
            [Pxx_trial, ~] = pwelch(cst1_detrend, h_win, noverlap, win_len, fsamp);
            [Pyy_trial, ~] = pwelch(cst2_detrend, h_win, noverlap, win_len, fsamp);

            % Cross-spectral accumulation
            Pxx_total = Pxx_total + Pxx_trial;
            Pyy_total = Pyy_total + Pyy_trial;
            Pxy_total = Pxy_total + Pxy_trial;
        end

        % Compute pooled multi-trial coherence after cross-trial accumulation
        % By accumulating the raw cross-spectra (Pxy_total) and
        % auto-spectra (Pxx_total, Pyy_total) across all trials
        % before applying the absolute square (abs(...).^2),
        % the cross-trial phase consistency is retained.
        % Formula: |Sum(Pxy)|^2 / (Sum(Pxx) * Sum(Pyy))
        denominator = Pxx_total .* Pyy_total;
        denominator(denominator == 0) = eps;

        msc_pooled = (abs(Pxy_total).^2) ./ denominator;
        perm_coherence(:, p) = msc_pooled;

    end %end of the "for p = 1:num_permutations" loop

    % Calculate true independent, disjoint data windows accumulated across all trials
    L_disjoint_per_trial = floor(num_samples / win_len);
    L_total_independent = L_disjoint_per_trial * num_trials; 

    % Bound individual permutations to avoid division-by-zero in atanh
    perm_coherence(perm_coherence >= 1) = 1 - eps; 

    % Transform each permutation to a Z-score first
    perm_z = sqrt(2 * L_total_independent) * atanh(sqrt(perm_coherence));

    % Take the average Z-scores across permutations
    coherence_z = mean(perm_z, 2);

    % Remove bias from Z-transformed coherence
    bias_idx = (f >= 250 & f <= 500);
    %raw_bias_level = mean(mean_coherence(bias_idx)); % remove raw bias
    z_bias_level = mean(coherence_z(bias_idx));

    % Rectify by subtracting the Z-scale noise floor
    coherence_rectified = coherence_z - z_bias_level;
    %coherence_rectified = mean_coherence - raw_bias_level; % raw coherence
    coherence_rectified(coherence_rectified < 0) = 0;

     %% Calculate raw confidence limit 
    alpha_95 = 0.05; % 95% Confidence level
    raw_threshold = 1 - (alpha_95^(1 / (L_total_independent - 1)));

    %% Transform the raw threshold into the Rosenberg Z-Score scale
    % Bound slightly below 1 to prevent atanh infinity errors
    if raw_threshold >= 1; raw_threshold = 1 - eps; end
    
    z_threshold = sqrt(2 * L_total_independent) * atanh(sqrt(raw_threshold));

    %% Rectify threshold against high-frequency noise floor
    rectified_threshold = z_threshold - z_bias_level;
    if rectified_threshold < 0; rectified_threshold = 0; end
    
    % Make CI_95_line the final significant line for 'coherence_rectified'
    CI_95_line = rectified_threshold; 

     %% Apply threshold gating for area under the curve calculation
    gated_coherence = coherence_rectified - CI_95_line;
    gated_coherence(gated_coherence < 0) = 0;

    %% Define frequency band indices
    delta_idx = (f >= delta_band(1) & f <= delta_band(2));
    alpha_idx = (f > 5   & f <= 15);
    beta_idx  = (f > 15  & f <= 35);
    gamma_idx = (f > 35  & f <= 60);

    %% Calculate threshold-gated AUC
    auc_delta = trapz(f(delta_idx), gated_coherence(delta_idx));
    auc_alpha = trapz(f(alpha_idx), gated_coherence(alpha_idx));
    auc_beta  = trapz(f(beta_idx),  gated_coherence(beta_idx));
    auc_gamma = trapz(f(gamma_idx), gated_coherence(gamma_idx));

    %color_idx = color_idx + 1;
end

% Plot coherence
figure;
plot(f, coherence_rectified, 'b', 'LineWidth', 1.5); hold on;
yline(CI_95_line, 'r--', '95% Confidence Limit', 'LineWidth', 1.5);
xlim([0 60]); % Limit view to delta, alpha, beta, and gamma bands
xlabel('Frequency (Hz)');
ylabel('Rectified Coherence (Z-Score Scale)');
title('Intramuscular Coherence with Significance Threshold');
grid on;

% Output the area under the curve for the largest MU pool size if looping
% through sub-pools
%auc_delta = auc_delta(end);
%auc_alpha = auc_alpha(end);
%auc_beta  = auc_beta(end);
%auc_gamma = auc_gamma(end);

end

% % --- VISUALIZATION BLOCK FOR THE LARGEST MU POOL SIZE ---
%     figure('Color', 'w', 'Name', sprintf('Condition %d Coherence Spectrum'));
%     hold on; grid on;
%
%     plot(f, coherence_rectified, 'b-', 'LineWidth', 2, ...
%         'DisplayName', sprintf('Max Pool (%d MUs Profile)', max_mu_per_pool));
%
%     L_per_trial = floor((num_samples - noverlap) / (win_len - noverlap));
%     L_total = L_per_trial * trials_counted;
%     conf_limit_msc = 1 - (0.01)^(1 / (L_total - 1));
%     conf_limit_rectified = max(0, conf_limit_msc - raw_bias_level);
%
%     yline(conf_limit_rectified, 'r--', 'Significance Threshold (\alpha=0.01)', ...
%         'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');
%
%     xlim([0 100]);
%     ylim([0 max(max(coherence_rectified)*1.1, conf_limit_rectified*1.5)]);
%     xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('Rectified Magnitude-Squared Coherence', 'FontSize', 12, 'FontWeight', 'bold');
%     title(sprintf('Pooled Intramuscular Coherence Spectrum (Condition %d)'), ...
%           'FontSize', 13, 'FontWeight', 'bold');
%
%     annotation_str = { ...
%         sprintf('\\bfDelta Area (0.5-5 Hz):\\rm %.3f', auc_delta(end)), ...
%         sprintf('\\bfAlpha Area (5-15 Hz):\\rm %.3f', auc_alpha(end)), ...
%         sprintf('\\bfBeta Area (15-35 Hz):\\rm %.3f', auc_beta(end)), ...
%         sprintf('\\bfGamma Area (35-60 Hz):\\rm %.3f', auc_gamma(end))};
%     text(55, max(ylim)*0.75, annotation_str, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10);
%
%     legend('Location', 'northeast');
%     drawnow;