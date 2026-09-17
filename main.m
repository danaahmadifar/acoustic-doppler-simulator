function main()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Akustik Projekt - Doppler-Effekt und Raumakustik-Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all;

% parameter
sim.fs             = 44100;
sim.frameSize      = 1024;
sim.nfft           = 2048;
sim.c              = 343;
sim.f0             = 220;
sim.gain           = 2.0;
sim.maxSpeed       = 25.0;
sim.raum           = [1000, 1000, 200];
sim.quellePos      = sim.raum / 2;
sim.kameraWinkel   = -pi / 4; % startwinkel der POV

sim.reflexionAktiv = true;
sim.daempfungAktiv = true;
sim.tonAktiv       = true;
sim.istPausiert    = false;

% hoerer startwerte
hoererPos = [150, 150, 40];
hoererVel = [0, 0, 0];
hoererYaw = pi / 4;

% ringpuffer
dsp.pufferLen    = round(10.0 * sim.fs);
dsp.puffer       = zeros(dsp.pufferLen, 1);
dsp.pufferPtr    = 1;
startDelay       = (norm(hoererPos - sim.quellePos) / sim.c) * sim.fs;
dsp.letztesDelay = repmat(startDelay, 7, 1);
dsp.letzterPan   = 0;

% audioquelle
isWav = false;
wavData = [];
wavPtr = 1;
synthPhase = 0;

win = hann(sim.frameSize);
fftFreqs = linspace(0, sim.fs / 2, sim.nfft / 2);
dt = sim.frameSize / sim.fs;
isRunning = true;

tasten = struct('w',0,'s',0,'a',0,'d',0,'space',0,'c',0,'q',0,'e',0);

% gui aufbauen
gui = gui_setup(sim);

set(gui.fig, 'KeyPressFcn',   @tasteRunter, 'KeyReleaseFcn', ...
    @tasteHoch, 'CloseRequestFcn', @beenden);
set(gui.chkRefl,   'Callback', @(s,~) setRefl(s.Value));
set(gui.chkDaempf, 'Callback', @(s,~) setDaempf(s.Value));
set(gui.sldSpeed,  'Callback', @(s,~) setSpeed(s.Value));
set(gui.sldGain,   'Callback', @(s,~) setGain(s.Value));
set(gui.sldFFT,    'Callback', @(s,~) setFFTMax(s.Value));
set(gui.btnPause,  'Callback', @(~,~) pauseWechsel());
set(gui.chkTon,    'Callback', @(s,~) setTon(s.Value));
set(gui.btnWav,    'Callback', @ladeWav);
set(gui.btnReset,  'Callback', @resetSinus);

% raumbox linien zeichnen
R = sim.raum;
boxX = [0 R(1) R(1) 0 0 0 R(1) R(1) 0 0 0 0 R(1) R(1) R(1) R(1)];
boxY = [0 0 R(2) R(2) 0 0 0 R(2) R(2) 0 0 R(2) R(2) 0 0 R(2)];
boxZ = [0 0 0 0 0 R(3) R(3) R(3) R(3) R(3) 0 0 0 0 R(3) R(3)];
set(gui.raumBox, 'XData', boxX, 'YData', boxY, 'ZData', boxZ);
xlim(gui.ax3D, [-20, R(1)+20]); 
ylim(gui.ax3D, [-20, R(2)+20]); 
zlim(gui.ax3D, [-10, R(3)+20]);

devWriter = audioDeviceWriter('SampleRate', sim.fs);
grafikZaehler = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulationsschleife
while isRunning && ishandle(gui.fig)

    % 1. hoerer bewegen
    cmdVor  = tasten.w - tasten.s;
    cmdDreh = tasten.a - tasten.d;
    cmdHoch = tasten.space - tasten.c;

    hoererYaw = hoererYaw + cmdDreh * 1.6 * dt;
    u_fwd     = [cos(hoererYaw), sin(hoererYaw), 0];

    zielVel   = u_fwd * (cmdVor * sim.maxSpeed) + [0, 0, cmdHoch * sim.maxSpeed * 0.4];
    hoererVel = 0.85 * hoererVel + 0.15 * zielVel;
    hoererPos = hoererPos + hoererVel * dt;

    hoererPos = max([5, 5, 2], min(sim.raum - 5, hoererPos));
    if hoererPos(3) <= 2.1 && hoererVel(3) < 0, hoererVel(3) = 0; end

    % kameradrehung mit Q und E
    sim.kameraWinkel = sim.kameraWinkel + (tasten.e - tasten.q) * 1.5 * dt;

    % 2. quellton generieren
    if sim.istPausiert || ~sim.tonAktiv
        quellBlock = zeros(sim.frameSize, 1);
    elseif isWav && ~isempty(wavData)
        idxEnd = wavPtr + sim.frameSize - 1;
        if idxEnd <= length(wavData)
            quellBlock = wavData(wavPtr:idxEnd);
            wavPtr = idxEnd + 1;
        else
            t1 = wavData(wavPtr:end);
            rest = sim.frameSize - length(t1);
            quellBlock = [t1; wavData(1:rest)];
            wavPtr = rest + 1;
        end
    else
        t = (0:sim.frameSize-1)' / sim.fs;
        quellBlock = sin(2 * pi * sim.f0 * t + synthPhase);
        synthPhase = mod(synthPhase + 2 * pi * sim.f0 * sim.frameSize / sim.fs, 2 * pi);
    end

    % 3. dsp abarbeiten
    [outBlock, metrics] = berechneDSP(quellBlock, hoererPos, hoererVel, u_fwd);
    devWriter(outBlock);

    % 4. grafik (~15 hz)
    grafikZaehler = grafikZaehler + 1;
    if grafikZaehler >= 3
        grafikZaehler = 0;

        % doppler verlauf
        gui.histShift = [metrics.shiftProzent; gui.histShift(1:end-1)];
        set(gui.linieShift, 'YData', gui.histShift, 'XData', 1:gui.numHist);

        % fft spektrum
        Y = abs(fft((outBlock(:, 1) + outBlock(:, 2)) * 0.5 .* win, sim.nfft));
        set(gui.linieFFT, 'XData', fftFreqs, 'YData', Y(1:sim.nfft / 2));

        % 3d objekte
        set(gui.hoerer, 'XData', hoererPos(1), 'YData', hoererPos(2), 'ZData', hoererPos(3));
        set(gui.linie,  'XData', [sim.quellePos(1), hoererPos(1)], ...
                        'YData', [sim.quellePos(2), hoererPos(2)], ...
                        'ZData', [sim.quellePos(3), hoererPos(3)]);

        mitte = (sim.quellePos + hoererPos) / 2;
        set(gui.txtDist, 'Position', [mitte(1), mitte(2), mitte(3) + 8], ...
                         'String', sprintf('%.1f m', metrics.dist));

        set(gui.pfeil, 'XData', hoererPos(1), 'YData', hoererPos(2), 'ZData', hoererPos(3), ...
                       'UData', u_fwd(1) * 60, 'VData', u_fwd(2) * 60, 'WData', 0);

        % kamera rotieren
        rMitte = sim.raum / 2;
        orbitR = norm(sim.raum(1:2)) * 0.95;
        kX = rMitte(1) + orbitR * cos(sim.kameraWinkel);
        kY = rMitte(2) + orbitR * sin(sim.kameraWinkel);
        kZ = sim.raum(3) * 1.5;
        set(gui.ax3D, 'CameraPosition', [kX, kY, kZ], 'CameraTarget', rMitte, ...
                      'CameraUpVector', [0, 0, 1], 'CameraViewAngle', 45);

        % hud text
        set(gui.hud, 'String', sprintf([ ...
            'Position: [%.0f, %.0f, %.0f] m\n', ...
            'Speed:    %.1f m/s\n', ...
            'Radial-v: %+.1f m/s\n', ...
            'Distanz:  %.1f m\n\n', ...
            'f_Send:   %.0f Hz\n', ...
            'f_Theo:   %.1f Hz\n', ...
            'Shift:    %+.2f %%\n\n', ...
            'Reflex:   %s\n', ...
            'Daempf:   %s'], ...
            hoererPos(1), hoererPos(2), hoererPos(3), ...
            norm(hoererVel), metrics.vRad, metrics.dist, ...
            sim.f0, metrics.fTheo, metrics.shiftProzent, ...
            textStatus(sim.reflexionAktiv), textStatus(sim.daempfungAktiv)));

        drawnow limitrate;
    end
end

release(devWriter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dsp berechnung
    function [outBlock, metrics] = berechneDSP(quellSignal, hPos, hVel, hFwd)
        bSize = length(quellSignal);
        sIndices = (1:bSize)';

        % in ringpuffer schreiben
        wI = mod((dsp.pufferPtr - 1 : dsp.pufferPtr + bSize - 2)', dsp.pufferLen) + 1;
        dsp.puffer(wI) = quellSignal * sim.gain;
        dsp.pufferPtr = mod(dsp.pufferPtr + bSize - 1, dsp.pufferLen) + 1;

        % quellen festlegen
        R_ = sim.raum; Q_ = sim.quellePos;
        if sim.reflexionAktiv
            nQ = 7;
            srcs = [ Q_; -Q_(1), Q_(2), Q_(3); 2*R_(1)-Q_(1), Q_(2), Q_(3); ...
                     Q_(1), -Q_(2), Q_(3); Q_(1), 2*R_(2)-Q_(2), Q_(3); ...
                     Q_(1), Q_(2), -Q_(3); Q_(1), Q_(2), 2*R_(3)-Q_(3) ];
        else
            nQ = 1;
            srcs = Q_;
        end

        mono = zeros(bSize, 1);
        for i = 1:nQ
            d = norm(hPos - srcs(i, :));
            tgtDly = max(15, (d / sim.c) * sim.fs);

            dlyRamp = dsp.letztesDelay(i) + (tgtDly - dsp.letztesDelay(i)) * ((sIndices - 1) / (bSize - 1));
            dsp.letztesDelay(i) = tgtDly;

            rPos = (dsp.pufferPtr - (bSize - sIndices + 1)) - dlyRamp;
            idx0 = floor(rPos);
            frac = rPos - idx0;

            s0 = dsp.puffer(mod(idx0 - 1, dsp.pufferLen) + 1);
            s1 = dsp.puffer(mod(idx0,     dsp.pufferLen) + 1);
            sample = s0 + frac .* (s1 - s0);

            if sim.daempfungAktiv, daempf = 25 / max(d, 25); else, daempf = 1; end
            if i == 1, mono = mono + sample * daempf;
            else,      mono = mono + sample * daempf * 0.25; end
        end

        % weiche peak-skalierung
        peakVal = max(abs(mono));
        if peakVal > 0.95
            mono = mono * (0.95 / peakVal);
        end

        % stufenloses stereo-panning
        relVec  = sim.quellePos - hPos;
        dist_   = norm(relVec);
        dir_    = relVec / max(dist_, 1e-4);
        u_right = [hFwd(2), -hFwd(1), 0];
        zielPan = max(-1, min(1, dot(u_right, dir_)));

        panRamp = dsp.letzterPan + (zielPan - dsp.letzterPan) * ((sIndices - 1) / (bSize - 1));
        dsp.letzterPan = zielPan;

        outBlock = [mono .* sqrt(0.5 * (1 - panRamp)), mono .* sqrt(0.5 * (1 + panRamp))];

        % metriken
        vRad_ = dot(hVel, dir_);
        metrics.vRad         = vRad_;
        metrics.dist         = dist_;
        metrics.fTheo        = sim.f0 * (1 + vRad_ / sim.c);
        metrics.shiftProzent = (vRad_ / sim.c) * 100;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gui anlegen
    function g = gui_setup(s)
        g.fig = figure('Name', '3D Akustik Labor', 'Position', [50, 50, 1300, 800], 'Color', [0.12 0.12 0.14]);

        % 3d raum
        g.ax3D = axes(g.fig, 'Position', [0.05, 0.40, 0.63, 0.54], 'Color', [0.05 0.05 0.07]);
        hold(g.ax3D, 'on'); grid(g.ax3D, 'on');
        set(g.ax3D, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
        xlabel(g.ax3D, 'X (m)', 'Color', 'w'); ylabel(g.ax3D, 'Y (m)', 'Color', 'w'); zlabel(g.ax3D, 'Z (m)', 'Color', 'w');

        g.raumBox = plot3(g.ax3D, nan, nan, nan, 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
        g.quelle  = plot3(g.ax3D, s.quellePos(1), s.quellePos(2), s.quellePos(3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        g.linie   = plot3(g.ax3D, nan, nan, nan, 'y--', 'LineWidth', 1.6);
        g.hoerer  = plot3(g.ax3D, nan, nan, nan, 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c');
        g.pfeil   = quiver3(g.ax3D, 0, 0, 0, 0, 0, 0, 0, 'Color', 'c', 'LineWidth', 2.2, 'MaxHeadSize', 1.5);
        g.txtDist = text(g.ax3D, 0, 0, 0, '', 'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold', ...
                         'HorizontalAlignment', 'center', 'BackgroundColor', [0.1 0.1 0.12]);

        % doppler verlauf
        g.axShift = axes(g.fig, 'Position', [0.05, 0.08, 0.30, 0.24], 'Color', [0.05 0.05 0.07]);
        hold(g.axShift, 'on'); grid(g.axShift, 'on');
        set(g.axShift, 'XColor', 'w', 'YColor', 'w');
        xlabel(g.axShift, 'Block', 'Color', 'w'); ylabel(g.axShift, 'Shift (%)', 'Color', 'w');
        title(g.axShift, 'Doppler-Shift \Delta f / f_0', 'Color', 'w');
        ylim(g.axShift, [-30, 30]); yline(g.axShift, 0, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

        g.numHist    = 80;
        g.histShift  = zeros(g.numHist, 1);
        g.linieShift = plot(g.axShift, g.histShift, 'c-', 'LineWidth', 2.0);

        % fft spektrum
        g.axFFT = axes(g.fig, 'Position', [0.38, 0.08, 0.30, 0.24], 'Color', [0.05 0.05 0.07]);
        hold(g.axFFT, 'on'); grid(g.axFFT, 'on');
        set(g.axFFT, 'XColor', 'w', 'YColor', 'w');
        xlabel(g.axFFT, 'Frequenz (Hz)', 'Color', 'w'); ylabel(g.axFFT, 'Betrag', 'Color', 'w');
        title(g.axFFT, 'FFT Spektrum', 'Color', 'w');
        xlim(g.axFFT, [0, 1500]); ylim(g.axFFT, [0, 80]);
        g.linieFFT = plot(g.axFFT, nan, nan, 'y-', 'LineWidth', 1.5);

        % steuerpanel
        g.panel = uipanel(g.fig, 'Title', 'Steuerung', 'Position', [0.70, 0.03, 0.28, 0.94], ...
                          'BackgroundColor', [0.15 0.15 0.18], 'ForegroundColor', 'w');

        mkCtrl = @(style, str, pos, val) uicontrol(g.panel, 'Style', style, 'String', str, ...
            'Position', pos, 'Value', val, 'BackgroundColor', [0.15 0.15 0.18], 'ForegroundColor', 'w');

        g.chkRefl   = mkCtrl('checkbox', 'Wandreflexionen',  [15, 715, 240, 20], s.reflexionAktiv);
        g.chkDaempf = mkCtrl('checkbox', 'Abstandsdaempfung', [15, 685, 240, 20], s.daempfungAktiv);

        mkCtrl('text', 'Max Speed (m/s):', [15, 645, 240, 16], 0);
        g.sldSpeed = uicontrol(g.panel, 'Style', 'slider', 'Min', 5, 'Max', 150, 'Value', s.maxSpeed, 'Position', [15, 625, 170, 20]);
        g.txtSpeed = mkCtrl('text', sprintf('%.0f m/s', s.maxSpeed), [190, 625, 65, 18], 0);

        mkCtrl('text', 'Lautstaerke (%):', [15, 585, 240, 16], 0);
        g.sldGain = uicontrol(g.panel, 'Style', 'slider', 'Min', 1, 'Max', 100, 'Value', 100, 'Position', [15, 565, 170, 20]);
        g.txtGain = mkCtrl('text', '100 %', [190, 565, 65, 18], 0);

        mkCtrl('text', 'FFT Max Pegel:', [15, 525, 240, 16], 0);
        g.sldFFT = uicontrol(g.panel, 'Style', 'slider', 'Min', 10, 'Max', 250, 'Value', 80, 'Position', [15, 505, 170, 20]);
        g.txtFFT = mkCtrl('text', '80', [190, 505, 65, 18], 0);

        g.btnPause = uicontrol(g.panel, 'Style', 'pushbutton', 'String', 'Pause [P]', 'Position', [15, 450, 240, 28], ...
                               'BackgroundColor', [0.3 0.3 0.35], 'ForegroundColor', 'w');
        g.chkTon   = mkCtrl('checkbox', 'Ton aktiv (An/Aus)', [15, 415, 240, 20], s.tonAktiv);
        set(g.chkTon, 'ForegroundColor', 'c');

        g.btnWav   = uicontrol(g.panel, 'Style', 'pushbutton', 'String', 'WAV laden...', 'Position', [15, 375, 115, 24]);
        g.btnReset = uicontrol(g.panel, 'Style', 'pushbutton', 'String', 'Sinus Reset',  'Position', [135, 375, 120, 24]);
        g.txtAudio = mkCtrl('text', 'Quelle: Sinus (220 Hz)', [15, 350, 240, 16], 0);
        set(g.txtAudio, 'ForegroundColor', [0.7 0.7 0.7]);

        g.hud = uicontrol(g.panel, 'Style', 'text', 'String', '', 'Position', [15, 15, 240, 315], ...
                          'BackgroundColor', [0.08 0.08 0.10], 'ForegroundColor', 'c', ...
                          'FontName', 'monospaced', 'FontSize', 9, 'HorizontalAlignment', 'left');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hilfsfunktionen
    function s = textStatus(v), if v, s = 'AN'; else, s = 'AUS'; end; end
    function setRefl(v),   sim.reflexionAktiv = (v == 1); end
    function setDaempf(v), sim.daempfungAktiv = (v == 1); end
    function setTon(v),    sim.tonAktiv = (v == 1); end
    function setSpeed(v),  sim.maxSpeed = v; set(gui.txtSpeed, 'String', sprintf('%.0f m/s', v)); end
    function setGain(v),   sim.gain = (v / 100) * 2.0; set(gui.txtGain, 'String', sprintf('%.0f %%', v)); end
    function setFFTMax(v), ylim(gui.axFFT, [0, v]); set(gui.txtFFT, 'String', sprintf('%.0f', v)); end
    function resetSinus(~,~), isWav = false; set(gui.txtAudio, 'String', 'Quelle: Sinus (220 Hz)'); end

    function pauseWechsel()
        sim.istPausiert = ~sim.istPausiert;
        if sim.istPausiert, set(gui.btnPause, 'String', 'Pausiert [P]', 'BackgroundColor', [0.7 0.2 0.2]);
        else,               set(gui.btnPause, 'String', 'Pause [P]',    'BackgroundColor', [0.3 0.3 0.35]); end
    end

    function ladeWav(~, ~)
        [datei, pfad] = uigetfile('*.wav', 'WAV Datei');
        if datei ~= 0
            [y, inFs] = audioread(fullfile(pfad, datei));
            if inFs ~= sim.fs, y = resample(y, sim.fs, inFs); end
            if size(y, 2) > 1, y = mean(y, 2); end
            wavData = y / (max(abs(y)) + eps);
            wavPtr = 1; isWav = true;
            set(gui.txtAudio, 'String', ['WAV: ', datei(1:min(12, length(datei)))]);
        end
    end

    function tasteRunter(~, e)
        switch lower(e.Key)
            case 'w', tasten.w = 1; case 's', tasten.s = 1;
            case 'a', tasten.a = 1; case 'd', tasten.d = 1;
            case 'space', tasten.space = 1; case 'c', tasten.c = 1;
            case 'q', tasten.q = 1; case 'e', tasten.e = 1;
            case 'p', pauseWechsel();
        end
    end

    function tasteHoch(~, e)
        switch lower(e.Key)
            case 'w', tasten.w = 0; case 's', tasten.s = 0;
            case 'a', tasten.a = 0; case 'd', tasten.d = 0;
            case 'space', tasten.space = 0; case 'c', tasten.c = 0;
            case 'q', tasten.q = 0; case 'e', tasten.e = 0;
        end
    end

    function beenden(~, ~), isRunning = false; delete(gui.fig); end
end