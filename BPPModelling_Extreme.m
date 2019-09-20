%% BPP Modelling extreme case
% Model to test stripped down case of very high order supergaussian beams to
% verify that conclusions drawn from only looking at n = 2--3.4 are valid.
format compact; clc
% Requires addaxis6 functions from "https://uk.mathworks.com/matlabcentral/fileexchange/9016-addaxis"

%% 1. Equations

SM_width_SG = @(x, curve, n)...                                            % Width of supergaussian via second moment
    sqrt(trapz(((x-0).^2).*curve, 2)/trapz(curve, 2)...                    % Second moment
    *(gamma(1/n)/gamma(3/n))*2^(2/n));                                     % Supergaussian adjustment == 4 for n == 2

%% 2. Orders and widths
% Fit curve to data to establish trend for supergaussian order extrapolation

% Load data
load DriverFWHMandP.mat
w0 = (1e-6)*FWHM/sqrt(2*log(2)); clear FWHM;                               % Beam waist exp(-2) [m]
n = p; clear p;                                                            % Supergaussian orders
w0 = w0(1:15);                                                             % Last five values had degraded shaping
n = n(1:15);
Num = 20;                                                                  % Number of beams to simulate

% Fit
[xData, yData] = prepareCurveData(w0/1e-6, n);
ft = fittype('a*(x^b) + 2', 'independent', 'x', 'dependent', 'y');         % Assume power law with intercept fixed at 2 (Gaussian)
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Lower = [0 1.5];                                                      % [a b c]
opts.Upper = [1e-5 5];
opts.StartPoint = (opts.Lower + opts.Upper)/2;
fitresult = fit(xData, yData, ft, opts);
nCoeffs = coeffvalues(fitresult);
nExtrapolation_equ = @(w0) (nCoeffs(1)*w0.^nCoeffs(2)) + 2;%nCoeffs(3);       % Make equation for calculating more n values

% Plot
F1 = figure(1); clf;
plot(fitresult, xData, yData);

% Define range of values to test subsequently
L = 3;                                                                     % Multiples of starting width to span over
w0_test = linspace(w0(1), L*w0(1), Num)/1e-6;                              % Beam widths to consider [um]
n_test = nExtrapolation_equ(w0_test);                                      % Corresponding n values
%n_test = 2*ones(size(w0_test));                                            % Override the test values

% Put extrapolated curve on the plot
figure(1); hold on;
plot(w0_test, n_test, '--', 'LineWidth', 1); hold off;
legend('Data', 'Fit', 'Extrapolated curve', 'Location', 'Best');
xlim([w0_test(1) w0_test(end)]); xlabel('Driver spot size [$\mu$m]');
ylim([1.8 max(n_test)]); ylabel('Supergaussian order');
F1.Children(2).Children(3).MarkerSize = 30;
F1.Children(2).Children(2).LineWidth = 3.5;
F1.Children(2).Children(1).LineWidth = 2.5;

%% 3. Constants

D = 0;                                                                     % Dipole status [1/0] for [on/off]
Plot = 0;                                                                  % Plots of source and detected profiles [1/0] for [on/off]
z = 15;                                                                    % Source to detector distance. Arbitrarily large to ensure far field [m]
lambda = 32e-9;                                                            % Radiation wavelength [m]
N = 2^12;                                                                  % Number of elements
x = -(N/2)+1:N/2;                                                          % Index axiss
dx = 1e-6;                                                                 % Source domain spatial resolution [m]
dX = (lambda*z)/(N*dx);                                                    % Detector domain spatial resolution [m]
alpha_s = 1.0367e-14; % -1.649e-14                                         % Short trajectory alpha for q = 25 [cm^2/W] at 4.7e14 Wcm^-2 [positive or negative?]
I_max = 4.7e14; %3.4e14                                                    % Peak focal intensity [W/cm^2]
qeff = 2.6;                                                                % Effective nonlinearity
Fresnel = ((w0_test*1e-6).^2)/(lambda*z);                                  % Fresnel number

%% 4. Define source term and propagate to far field

SG = zeros(length(w0_test), N);                                            % Supergaussian driver foci
SGh = zeros(length(w0_test), N);                                           % Supergaussian harmonic sources
SGff = zeros(length(w0_test), N);                                          % Detected harmonic profiles
wq = zeros(length(w0_test), 1);                                            % Harmonic source size
Wz = zeros(length(w0_test), 1);                                            % Detected harmonic profile size

L1 = 1e-4;                                                                 % Limits for plots if they're used
L2 = 10e-3;
for i = 1:length(w0_test)
    k = (Num-i)/Num;                                                       % Colour index
    
    % Source - IR
    SG(i,:) = exp(-abs(x*dx/(w0_test(i)*dx)).^n_test(i));                  % Supergaussian FIELD demagnified into secondary source plane
    % Source - Harmonic
    SGh(i,:) = (SG(i,:).^qeff).*exp(-1i*D*alpha_s*I_max*abs(SG(i,:)).^2);  % Harmonic FIELD
    wq(i) = SM_width_SG(x*dx, abs(SGh(i,:)).^2, n_test(i))/1e-6;           % Harmonic secondary source width [um]

    % Propagate harmonic sources to far field
    SGff(i,:) = fftshift(fft(fftshift(SGh(i,:))))*dx;                      % Propagated field
    Wz(i) = SM_width_SG(x*dX, abs(SGff(i,:)).^2, 2)/1e-6;                  % Harmonic secondary source width [um]

    if Plot == 1
        % Plot source terms
        F2 = figure(2); hold on;
        title('Harmonic source');
        yyaxis left
        plot(x*dx, abs(SGh(i,:)).^2, '-', 'Color', [k k k]);               % Harmonic source intensities
        yyaxis right
        plot(x*dx, unwrap(angle(SGh(i,:))), '-');                          % Harmonic source phases
        xlim([-L1 L1]);

        % Plot detected profiles
        F3 = figure(3); hold on;
        title('Far field harmonic profiles');
        yyaxis left
        plot(x*dX, abs(SGff(i,:)).^2, '-', 'Color', [k k k]);              % Detected harmonic intensities
        yyaxis right
        plot(x*dX, unwrap(angle(SGff(i,:))), '-');                         % detected harmonic phases
        xlim([-L2 L2]);
        
        F2.Children.YAxis(1).Color = [0 0 0];
        F3.Children.YAxis(1).Color = [0 0 0];

    else
    end
end

BPP = wq.*(Wz/(z*1e6));                                                    % Calculate the BPP once simulations complete


%% 5. Plot BPP, wq and Wz

if n_test(2) == 2
    set(0, 'defaultaxesfontsize', 35);
    F4 = figure(4); clf;
    yyaxis left;
    plot(w0_test, BPP/1e-3, 'Color', [0 0 0.8]); ylim([0.5 1.5]*BPP(1)/1e-3); ylabel('BPP [$\mu$m$\:$mrad]');
    addaxis(w0_test, Wz/1e3, '-', 'Color', [0 0.6 0]); addaxislabel(2, 'W(z) [mm]');
    yyaxis right;
    plot(w0_test, wq, 'Color', [0.8 0 0]);  ylabel('$w_q$ [$\mu$m]');
    xlim([w0_test(1) w0_test(end)]); xlabel('$w_0$ [$\mu$m]');

    % Color the axes and lines
    F4.Children(1).YAxis(1).Color = [0 0 0.8];                             % BPP
    F4.Children(1).YAxis(2).Color = [0.8 0 0];                             % wq
    F4.Children(2).YAxis.Color = [0 0.6 0];                                % Wz

    % Space the axes nicely
    Pos = F4.Children(1).Position;
    F4.Children(1).Position = Pos - [0 -0.03 0.1 0];
    F4.Children(2).Position = Pos + [0 0.03 0.02 0];                       % [left edge, bottom edge, right edge, top edge]
    F4.Children(2).XAxis.Visible = 'off';

else                                                                       % If n varies...
    
    set(0, 'defaultaxesfontsize', 35);
    F5 = figure(5); clf;
    yyaxis left;
    plot(w0_test, BPP/1e-3, 'Color', [0 0 0.8]); ylim([0.5 1.5]*BPP(1)/1e-3); ylabel('BPP [$\mu$m$\:$mrad]');
    addaxis(w0_test, Wz/1e3, '-', 'Color', [0 0.6 0]); addaxislabel(2, 'W(z) [mm]');
    yyaxis right;
    plot(w0_test, wq, 'Color', [0.8 0 0]);  ylabel('$w_q$ [$\mu$m]');
    xlim([w0_test(1) w0_test(end)]); xlabel('$w_0$ [$\mu$m]');

    % Color the axes and lines
    F5.Children(1).YAxis(1).Color = [0 0 0.8];                             % BPP
    F5.Children(1).YAxis(2).Color = [0.8 0 0];                             % wq
    F5.Children(2).YAxis.Color = [0 0.6 0];                                % Wz

    % Space the axes nicely                                                % [left edge, bottom edge, right edge, top edge]
    Pos = F5.Children(1).Position;
    F5.Children(1).Position = Pos - [0 -0.03 0.1 0];
    F5.Children(2).XAxis.Visible = 'off';
    Pos_new = F5.Children(1).Position;
    
    % Make new axes to fit n onto upper x axis
    ax2 = axes('Position', Pos_new, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
    ax2.YAxis.Visible = 'off';
    F5.Children(1).XAxis.Limits = [n_test(1) n_test(end)];

    % More axis adjustments to get everything to fit into your right hand monitor
    F5.Children(1).XAxis.TickValues = linspace(n_test(1), n_test(end), length(F5.Children(2).XAxis.TickValues));
    F5.Children(1).XAxis.TickLabelFormat = '%.1f';
    F5.Children(1).XLabel.String = '$n$';
    F5.Children(2).Position = Pos_new - [0 0 0 0.08];                     % [left edge, bottom edge, right edge, top edge]
    F5.Children(1).Position = Pos_new + [0 0 0 -0.08];
    F5.Children(3).Position = Pos_new - [0 0 -0.12 0.08];
    xlim([n_test(1) n_test(end)]);

end

%% 6. Plot variation of the supergaussian adjustment (text line 3, section 1)

Adjustment = @(n) 2.^(2./n).*gamma(1./n)./gamma(3./n);                     % Factor in equation to determine supergaussian beam size
WidthDenominator = @(n, nl) nl.^(1./n);                                    % Denominator in equation to calculate harmonic source width

F6 = figure(6); clf;
yyaxis left
plot(n_test, Adjustment(n_test), 'Color', [0 0 0.8]);
ylabel('Second moment adjustment');
yyaxis right
plot(n_test, WidthDenominator(n_test, qeff), 'Color', [0.8 0 0]);
ylabel('Dipole model width denominator');

xlabel('$n$');
F6.Children.YAxis(1).Color = [0 0 0.8];
F6.Children.YAxis(2).Color = [0.8 0 0];
axis tight;

% Note, Dipole model source size denominator (red curve) tends toward unity
% as n tends to \infty, so the faster-than-linear increase in wq with Wz
% exists only at lowish n values. This makes sense because the difference
% in source profile between n=2--5 is a much larger change than from
% n=10--15

% Red curve essentially describes how as the driver profile becomes a
% uniform intensity tophat, the source will eventually match its size. Its
% asymptotic approach to =1 is can be seen in figure 5 (with n allowed to
% vary) where the red curve for wq has a slight bend in the low n region.

F7 = figure(7); clf;
yyaxis left
plot(n_test, WidthDenominator(n_test, qeff), 'Color', [0.8 0 0]);
ylabel('$\sqrt[n]{q_\mathrm{eff}}$');
yyaxis right
plot(n_test(1:end-1), diff(WidthDenominator(n_test, qeff)), 'Color', [0 0 0.8]);
ylabel('$\frac{\mathrm{d}}{\mathrm{d}n}\sqrt[n]{q_\mathrm{eff}}$');
xlabel('$n$');

F7.Children.YAxis(1).Color = [0.8 0 0];
F7.Children.YAxis(2).Color = [0 0 0.8];
axis tight;

%% 7. Plot rates of decrease of Wz
% This won't plot until you've run this script 4 times with D = [0 1] and 
% n = [varies fixed]

set(0, 'defaultaxesfontsize', 40);
try
    if (D == 1) && (n_test(2) == 2)
        Wz_D1_nF = Wz;
    else
        if (D == 1)
            Wz_D1_nV = Wz;
        else
            if (n_test(2) == 2)
                Wz_D0_nF = Wz;
            else
                Wz_D0_nV = Wz;
            end
        end
    end
    
    F8 = figure(8);
    plot(w0_test, Wz_D1_nF./Wz_D0_nF, 'Color', [0 0 0.8]); hold on;
    plot(w0_test, Wz_D1_nV./Wz_D0_nV, 'Color', [0.8 0 0]); hold off;
    xlabel('$w_0$ [$\mu$m]');
    ylabel('$\frac{W(z)\mbox{ dipole on}}{W(z)\mbox{ dipole off}}$');
    legend('$n$ fixed', '$n$ varies', 'Location', 'Best');
    xlim([w0_test(1) w0_test(end)]);
    ylim([1.25*max(Wz_D1_nF./Wz_D0_nF) 0.75*min(Wz_D1_nV./Wz_D0_nV)]);
    
catch ME

end
    
