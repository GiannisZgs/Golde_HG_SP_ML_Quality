function [H_w,H_s,ZCR1,ZCR2,SF,Compr] = waveform_quality_scores(beats)
    % Calculates multiple ECG signal quality metrics
    % Returns waveform entropy, spectral entropy, zero-crossing rates,
    % spectral flatness, and compressibility measures
    
    H_w = waveform_entropy(beats, 100);
    H_s = spectral_entropy(beats, 100);
    ZCR1 = zero_crossing_rate(beats);
    ZCR2 = num_local_extrema(beats);
    SF = spectral_flatness(beats);
    Compr = compressed_size(beats);
    
end

function H = waveform_entropy(beats, nbins)
    [~, N] = size(beats);
    H = zeros(1, N);
    for i = 1:N
        beat = beats(:,i);
        [counts, ~] = histcounts(beat, nbins, 'Normalization', 'probability');
        counts = counts + 1e-10;
        H(i) = -sum(counts .* log2(counts));
    end
end

function H = spectral_entropy(beats, nbins)
    X = abs(fft(beats,64)); %operate along columns
    X = X(2:floor(end/2),:);  % remove DC and mirror
    [~, N] = size(beats);
    H = zeros(1, N);
    for i = 1:N
        spec = X(:,i);
        [counts, ~] = histcounts(spec, nbins, 'Normalization', 'probability');
        counts = counts + 1e-10;
        H(i) = -sum(counts .* log2(counts));
    end
end

function Z = zero_crossing_rate(beats)
    %1st order ZC
    [~, N] = size(beats);
    Z = zeros(1, N);
    for i = 1:N
        b = beats(:,i);
        zc = diff(sign(b));
        Z(i) = sum(zc ~= 0);
    end
end

function SF = spectral_flatness(beats)
    [~, N] = size(beats);
    SF = zeros(1, N);
    for i = 1:N
        b = beats(:,i);
        X = abs(fft(b));
        X = X(2:floor(end/2));  % remove DC and mirror
        geo_mean = exp(mean(log(X + 1e-10)));
        arith_mean = mean(X);
        SF(i) = geo_mean / arith_mean;
    end
end

function E = num_local_extrema(beats)
    %2nd order derivative
    [~, N] = size(beats);
    E = zeros(1, N);
    for i = 1:N
        b = beats(:,i);
        d = diff(b);
        E(i) = sum(diff(sign(d)) ~= 0);
    end
end

function C = compressed_size(beats)
    [~, N] = size(beats);
    C = zeros(1, N);
    for i = 1:N
        beat = beats(:,i);
        temp = [tempname '.mat'];
        save(temp, 'beat', '-mat');
        fileInfo = dir(temp);
        C(i) = fileInfo.bytes;
        delete(temp);
    end
end

