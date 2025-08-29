function features = extract_ecg_heartbeat_features(beat, fs)

    beat = beat(:);
    L = length(beat);
    
    % --- Peak-to-peak amplitude ---
    features.peak_to_peak = max(beat) - min(beat);

    % --- R-peak detection (largest absolute peak) ---
    [~, R_idx] = max(abs(beat));

    % --- Q and S detection (local extrema around R) ---
    
        
    search_range = round(0.1 * fs); % ±100 ms
    q_range = max(1, R_idx-search_range):R_idx;
    s_range = R_idx:min(L, R_idx+search_range);
    if beat(R_idx) < 0
        %Inverted R-peak
        [~, q_rel_idx] = max(beat(q_range));
        [~, s_rel_idx] = max(beat(s_range));
    else
        [~, q_rel_idx] = min(beat(q_range));
        [~, s_rel_idx] = min(beat(s_range));
    end
    q_idx = q_range(q_rel_idx);
    q_amp = beat(q_idx);
    if q_amp > 0
        q_amp = -q_amp;
    end
    s_idx = s_range(s_rel_idx);

    if ~isempty(q_idx) && ~isempty(s_idx)
        features.qrs_duration = (s_idx - q_idx) / fs;
    else
        features.qrs_duration = NaN;
    end

    % --- P-wave detection (local extremum before Q) ---
    p_window =1:q_idx; % max(1, q_idx - round(0.2 * fs))
    if beat(R_idx) < 0
        %Inverted R-peak
        [~, p_rel_idx] = min(beat(p_window));
    else
        [~, p_rel_idx] = max(beat(p_window));
    end
    p_idx = p_window(p_rel_idx);
    features.p_wave_amp = abs(beat(p_idx) - q_amp);

    % --- T-wave detection (local extremum after S) ---
    t_window = s_idx:L;
    if beat(s_idx) < 0
        [~, t_rel_idx] = max(beat(t_window));
    else
        [~, t_rel_idx] = min(beat(t_window));
    end
    t_idx = t_window(t_rel_idx);
    features.t_wave_amp = abs(beat(t_idx) - q_amp);
    
    % Optionally: return indices as well
    features.indices = struct('P', p_idx, 'Q', q_idx, 'R', R_idx, 'S', s_idx, 'T', t_idx);
end
