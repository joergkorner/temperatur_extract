# temperatur_extrakt.py — Temperaturwerte aus IGC-Logs, ohne Spur

Was es tut: Aus jedem Flug-Log werden je Zeitfenster und 25-hPa-Schicht eine Temperatur (aus Druckhöhe gegen GPS-Höhe,
hypsometrische Gleichung) und ein paar Begleitzahlen gezogen. Keine Spur, keine Positionen feiner als 1 km, keine
Dateinamen, kein Pilot. Methode: Korner 2026, https://doi.org/10.31223/X5PJ7S.

Lauf (neben dem Segmentierer, auf demselben Archiv; nur numpy nötig):

    python3 temperatur_extrakt.py IGC_ORDNER AUSGABE/alpen --kerne 8

Testlauf vorher:

    python3 temperatur_extrakt.py IGC_ORDNER AUSGABE/test --probe 1000

Ergebnis: zwei gepackte Tabellen, AUSGABE/alpen_fluege.csv.gz (eine Zeile je Log) und AUSGABE/alpen_schichten.csv.gz
(eine Zeile je Log, Fenster und Schicht). Die Spalten stehen oben im Skript. Grössenordnung: 950'000 Logs, Laufzeit rund
eine Stunde auf 8 Kernen, Ausgabe einige hundert MB gepackt. Unterordner werden mitgenommen.

Kennung: Standard ist --id hash (10 Zeichen aus dem Dateinamen, nicht zurückrechenbar). Dieselbe Kennung wie beim
Segmentierer wäre ideal, falls Du dort --id mit einem eigenen Schlüssel benutzt; sonst passen die Zeilen über Datum,
Sekunde des ersten Fixes und Position des ersten Fixes zusammen (t_first, lat_first, lon_first).

Kaputte Dateien werden übersprungen und auf stderr gemeldet; Logs mit nur einer Höhe bekommen nur die Zeile in
_fluege (two_altitudes = 0), damit das Archiv vollständig beschrieben ist.
