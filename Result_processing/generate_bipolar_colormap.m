function cmap = generate_bipolar_colormap(lutsize, neutral, interpMethod)
    %lutsize = 256;
    %neutral = 0; 
    %interpMethod = 'cubic';

    % Generates a bipolar colormap

    if neutral < 0.5
        data = [
            0, 1, 1; % cyan
            0, 0, 1; % blue
            neutral, neutral, neutral; % neutral
            1, 0, 0; % red
            1, 1, 0; % yellow
        ];
    else
        data = [
            0, 0, 1; % blue
            0, 1, 1; % cyan
            neutral, neutral, neutral; % neutral
            1, 1, 0; % yellow
            1, 0, 0; % red
        ];
    end

    xi = linspace(0, 1, size(data, 1));
    xnew = linspace(0, 1, lutsize);

    cmap = interp1(xi, data, xnew, interpMethod);
    cmap = max(0, min(cmap, 1)); % Ensure values are within [0, 1]
    cmap = flip(cmap);
end