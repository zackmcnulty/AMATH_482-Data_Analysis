% gives plots of the common windows used in Gabor Transforms
tspan = 0:0.01:10;
tau = max(tspan) / 2;
a = 0.2;

% gaussian
wg = exp(-a*(tspan - tau).^2);

% super gaussian
wsg = exp(-a*(tspan - tau).^10);

% Shannon
ws = abs(tspan - tau) <= 1/(2*a);

% Mexican Hat
wmh = (1-(tspan - tau).^2).*exp(-(tspan - tau).^2 ./ 2);
%wmh = (2 ./ (sqrt(3*a)*pi^0.25)) .* (1 - ((tspan - tau) ./ a).^2).*exp(-1.*(tspan - tau).^2 ./ (2*a^2));
    

figure (1)
plot(tspan, [wg.' wsg.' ws.' wmh.'], 'linewidth', 5);
hleg = legend({"Gaussian", "Super Gaussian", "Shannon", "Mexican Hat"});
title(hleg, "Window Types")
set(gca, 'fontsize' ,25);
xlabel("Time")
ylabel("Amplitude")
title("Comparison of Different Gabor Windows")
ylim([-1, 2])


%% Plot the spectrogram example



tspan = 0:0.01:pi;
y = sin(5*tspan.^2);

figure(1)
spectrogram(y, 'yaxis');
ylim([0,0.1])

figure(2)
plot(tspan, y);
ylim([-1.5, 1.5])







