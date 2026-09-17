# Real-Time Acoustic & Doppler Simulator

Echtzeit-Simulation von akustischen Ausbreitungseffekten, Doppler-Frequenzverschiebungen und räumlichem Stereo-Hören in MATLAB.
Es können eigene .wav-Dateien genutzt werden, um die Effekte mit dieser zu testen.

## Features
* **Echtzeit-Doppler-Effekt:** Kontinuierliche Abtastung über zeitvariable Delay-Rampen und lineare Puffer-Interpolation.
* **3D-Stereo-Panning:** Automatische Pegelaufteilung auf linken und rechten Kanal anhand des Kopfausrichtungsvektors.
* **Wandreflexionen:** Simulation von 6 Bildquellen an den Raumgrenzen.
* **Live-Analyse:** Parallele Visualisierung von Flugpfad, Distanz, theoretischem Doppler-Shift (%) und FFT-Ausgangsspektrum.

## Steuerung des "Helicopterss"
* W / S: Vorwärts / Rückwärts
* A / D: Drehen
* Space / C: Steigen / Sinken
* Q / E: 3D-Kamera rotieren
* P: Audio-Quelle pausieren

## Voraussetzungen
* MATLAB (getestet auf R2026b)
* Audio Toolbox ("audioDeviceWriter")

## Starten
Einfach "main.m" in MATLAB öffnen und ausführen:# acoustic-doppler-simulator
Ein kleines Projekt, mit dem man den Doppler Effekt und und Wandreflexionen in einem Raum simulieren kann.
