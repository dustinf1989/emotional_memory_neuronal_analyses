function [width, ratio] = get_spike_width_ratio(mean_wav)
% Get spike width and spike ratio from an average, unormalized waveform
    [mn, imin] = min(mean_wav);
    [mx, imax] = max(mean_wav);
    width = abs(imax - imin) / length(mean_wav);
    if width > 0.75
        disp('Check spike width')
        [mn, imin] = min(mean_wav(5:end-5));
        [mx, imax] = max(mean_wav(5:end-5));
        imin = imin + 5;
        imax = imax + 5;
        width = abs(imax - imin) / length(mean_wav);
    end

%     peak_valley = abs(mx) + abs(mn);

%     peak = mean_wave[maxi]
%     valley = mean_wave[mini] % wc in Leibold code

    ratio = abs(mx) / abs(mn);
end