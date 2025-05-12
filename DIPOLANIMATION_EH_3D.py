import numpy as np
import pyvista as pv
import math
import os
import shutil
import sys
import imageio.v2 as imageio
import re


if os.path.exists("DipolAnimation"):
    shutil.rmtree("DipolAnimation")  # Vorheriges Verzeichnis löschen, falls vorhanden
os.mkdir("DipolAnimation")  # Verzeichnis für Animation erstellen

save_gif = True  # Bei Bedarf auf True setzen, um GIF zu speichern


# Physikalische Konstanten und Parameter
_c = 299792458  # Lichtgeschwindigkeit
_c2 = 8.9876e16  # Lichtgeschwindigkeit^2
_epsilon0 = 8.8542e-12
_my0 = 1.2566e-6
_Z0 = 376.7303134

# Dipol-Parameter
_p0 = 100.0  # Amplitude des Dipolvektors
_Wellenlaenge = 256.0  # Wellenlänge (Simulationseinheiten)
_Lamda_viertel = _Wellenlaenge / 4.0
_w = 2 * np.pi * _c / _Wellenlaenge  # Kreisfrequenz
_Periode = 100  # Anzahl Zeitschritte pro Periode
_animation_duration = 3.0  # Animationsdauer in Sekunden für eine Periode
_T = _Wellenlaenge / _c  # Periodendauer
_dt = _T / _Periode  # Zeitschritt

# Darstellungs-Flags
Hertzdipol = True
Stabdipol = False
E_Feldlinien_Flag = True  # elektrische Feldlinien anzeigen
H_Feldlinien_Flag = True  # magnetische Feldlinien (H-Feld) anzeigen
Energiestrom_Flag = False  # Energiestrom-Pfeile anzeigen
EnergieMakro_Flag = False  # Makro-Energiestromdichte (wenige, große Pfeile)

# Skalierungsfaktoren
_Groessenfaktor = 3
_Horizontalfaktor = 2.0 / _Groessenfaktor
_Vertikalfaktor = 2.0 / _Groessenfaktor

# Parameter für Feldlinien-Berechnung
_Punktzahl_max = 3500
_Grobfaktor2 = 0.25  # Schrittweitenfaktor für Feldlinienintegration
_ZahlDerRechenschritte = int(_Punktzahl_max / _Grobfaktor2 + 100)
_Linienabst = 8.0  # Abstand zwischen Feldlinien
_Linien_pro_Welle = _Wellenlaenge / 4.0 / _Linienabst

# Parameter für magnetische Feldlinien (H-Feld) in der xy-Ebene
_radd = [
    17,
    5.0,
    16.0,
    26.0,
    35.0,
    43.0,
    50.0,
    55.0,
    60.0,
    64.0,
    68.0,
    73.0,
    79.0,
    86.0,
    94.0,
    103.0,
    112.0,
    125.0,
]

# Parameter für Energiestrom-Pfeile (Abstände und Pfeilgröße)
if EnergieMakro_Flag:
    _xAbst = 40.0  # Pfeil-Abstand in Pixel (Original)
    _yAbst = 40.0
    _PfeilPktzahl = 20  # Pfeilgrößenfaktor (20 für Makro)
else:
    _xAbst = 20.0
    _yAbst = 20.0
    _PfeilPktzahl = 10  # Pfeilgrößenfaktor (10 für Mikro)

# Farbgrenzen und Farben für Poynting-Vektor (Energiestromdichte)
mmm = 5.0
Farbgrenze0 = mmm * 1200 * 1e-10
Farbgrenze1 = mmm * 600 * 1e-10
Farbgrenze2 = mmm * 300 * 1e-10
Farbgrenze3 = mmm * 150 * 1e-10
Farbgrenze4 = mmm * 75 * 1e-10
Farbgrenze5 = mmm * 10 * 1e-10

Farbe_S0 = "#FF0000"
Farbe_S1 = "#CC0000"
Farbe_S2 = "#AA1111"
Farbe_S3 = "#881188"
Farbe_S4 = "#555555"
Farbe_S5 = "#999999"

# Farben für Feldlinien
Farbe_fieldline_pos = "#008800"  # 008800"  # grün
Farbe_fieldline_neg = "#0000AA"  # blau

print(
    "Important: Currently this application creates a directory named 'DipolAnimation' which contains the images for each timestep. They need to be converted into a gif afterwards"
)


# E-Feld berechnen
def E_berechnen(t, x, y):
    if Hertzdipol:
        r2 = x * x + y * y
        r = np.sqrt(r2)
        if r == 0:
            return (0.0, 0.0, 0.0)
        tt = t - r / _c
        p = _p0 * np.cos(_w * tt)
        pd1 = -_w * _p0 * np.sin(_w * tt)
        cosa = y / r
        Efy = (_w**2 * p) / (_c2 * r)
        Ep = Efy - pd1 / (_c * r2) - p / (r2 * r)
        Er = 3 * pd1 / (_c * r2) - Efy + 3 * p / (r2 * r)
        Er *= cosa
        Ex = Er * (x / r)
        Ey = Er * cosa + Ep
        E = np.sqrt(Ex * Ex + Ey * Ey)
        return (Ex, Ey, E)
    elif Stabdipol:
        Lplus = y + _Lamda_viertel
        Lminus = y - _Lamda_viertel
        rplus = np.sqrt(x * x + Lplus * Lplus)
        rminus = np.sqrt(x * x + Lminus * Lminus)
        if rplus == 0 or rminus == 0:
            return (0.0, 0.0, 0.0)
        tt1 = t - rplus / _c
        tt2 = t - rminus / _c
        sinplus = (1.0 / rplus) * np.cos(_w * tt1)
        sinminus = (1.0 / rminus) * np.cos(_w * tt2)
        Ey = (sinplus + sinminus) / 40.0  # _Stabfaktor = 40
        if x == 0:
            Ex = 0.0
        else:
            Ex = -((Lplus / x * sinplus) + (Lminus / x * sinminus)) / 40.0
        E = np.sqrt(Ex * Ex + Ey * Ey)
        return (Ex, Ey, E)


def H_berechnen(t, x, y):
    """Berechnet das magnetische Feld (H-Phi-Komponente) an Punkt (x,y) zur Zeit t."""
    H = 0.0
    if Hertzdipol:
        r2 = x * x + y * y
        r = np.sqrt(r2)
        if r == 0:
            return 0.0
        tt = t - r / _c
        p = _p0 * np.cos(_w * tt)
        pd1 = -_w * _p0 * np.sin(_w * tt)
        H = pd1 / (_Z0 * r) + (_w * p) / (_Z0 * _c * r)
    elif Stabdipol:
        # Vereinfachte Berechnung für Stabdipol (nur Betrag, H in xy-Ebene)
        r = np.sqrt(x * x + y * y)
        if r == 0:
            return 0.0
        # angenommen H proportional zu current_I und ~1/r
        current_I = _p0 * np.sin(_w * t)  # Strom im Dipol (vereinfachtes Modell)
        H = current_I / (2 * np.pi * r)
    return H


def S_berechnen(t, x, y):
    """Berechnet den Poynting-Vektor S=(Sx, Sy) und |S| in der xy-Ebene an Punkt (x,y) zur Zeit t."""
    # In dieser Simulation: S = E x H (mit E- und H-Feldkomponenten).
    Ex, Ey, _ = E_berechnen(t, x, y)
    H_phi = H_berechnen(t, x, y)
    # E-Feld liegt in r-Ebene (Ex, Ey), H-Feld steht senkrecht dazu (aus der Ebene heraus für radiales E)
    # Poynting zeigt radial vom Dipol weg: S = E x H * ê_richtung
    # Vereinfachung: Sx = Ey * H_phi (entspricht z-Komponente von E x H), Sy = -Ex * H_phi.
    Sx = Ey * H_phi  # positiv nach außen (x-Richtung)
    Sy = -Ex * H_phi  # positiv nach außen (y-Richtung)
    S = np.sqrt(Sx * Sx + Sy * Sy)
    return (Sx, Sy, S)


# Hilfsfunktionen für Feldlinien-Berechnung
def F(x_val, t):
    """Hilfsfunktion F(x) für Hertzdipol (Nullstellenbestimmung auf x-Achse)."""
    zw = -_w * (x_val / _c - t)
    return np.cos(zw) - (x_val * _w / _c) * np.sin(zw)


def FSt(x_val, t):
    """Hilfsfunktion FSt(x) für Stabdipol."""
    rpm2 = np.sqrt(x_val * x_val + _Lamda_viertel * _Lamda_viertel)
    if rpm2 == 0:
        return 0.0
    tt = t - rpm2 / _c
    return (2 / rpm2) * np.cos(_w * tt)


def Nullstelle(x1, x2, t, func):
    """Bisektionsverfahren zur Nullstellensuche von func im Intervall [x1, x2]."""
    tol = 0.001
    max_iter = 25
    f1 = func(x1, t)
    f2 = func(x2, t)
    x_mid = None
    if f1 == 0:
        return x1
    if f2 == 0:
        return x2
    for _ in range(max_iter):
        x_mid = (x1 + x2) / 2.0
        f_mid = func(x_mid, t)
        if abs(x2 - x1) <= tol or f_mid == 0:
            break
        if f2 * f_mid > 0:
            # Nullstelle liegt in [x1, x_mid]
            x2 = x_mid
            f2 = f_mid
        else:
            # Nullstelle liegt in [x_mid, x2]
            x1 = x_mid
            f1 = f_mid
    return x_mid if x_mid is not None else (x1 + x2) / 2.0


# Bestimme "Grenzlinien" (Nullstellen) auf der x-Achse, die Feldlinien trennen
Zahl_der_Extrema = 0
_Grenzlinie = []


def Grenzlinien_bestimmen(t):
    global Zahl_der_Extrema, _Grenzlinie
    _Grenzlinie = []
    if Hertzdipol:
        Keulenzahl = 7
        _Groessenfaktor2_val = 2 if False else _Groessenfaktor  # Pot3D_Flag nicht aktiv
        Schrittweite = _Lamda_viertel / 4.0
        x1 = _Linienabst / 2.0 + 3.0  # Start etwas jenseits Nahgrenze (Nahgrenze = 3)
        x2 = x1 + Schrittweite
        ij = 0
        Fensterende = False
        while not Fensterende:
            f1 = F(x1, t)
            f2 = F(x2, t)
            if f1 * f2 <= 0:
                nullst = Nullstelle(x1, x2, t, F)
                _Grenzlinie.append(nullst)
                ij += 1
                Zahl_der_Extrema = ij
            x1 = x2
            x2 = x2 + Schrittweite
            if ij > 0:
                if (
                    _Grenzlinie[ij - 1]
                    > _Groessenfaktor2_val * Keulenzahl * _Lamda_viertel
                ):
                    Fensterende = True
            if x1 > _Groessenfaktor2_val * Keulenzahl * _Lamda_viertel:
                Fensterende = True
            if ij >= 100:  # Sicherheit
                Fensterende = True
    if Stabdipol:
        _Grenzlinie = []
        ij = 0
        tt = t
        # Keulenzahl abhängig vom Zoomfaktor
        if _Groessenfaktor < 1:
            Keulenzahl = 4
        elif _Groessenfaktor == 1:
            Keulenzahl = 5
        elif _Groessenfaktor == 2:
            Keulenzahl = 7
        elif _Groessenfaktor == 4:
            Keulenzahl = 12
        else:
            Keulenzahl = 7
        for n in range(Keulenzahl):
            # Normiere Zeitphase in [0, π/_w]
            while tt > np.pi / _w:
                tt -= np.pi / _w
            aa1 = (2 * n) * _Lamda_viertel + _c * tt
            # (Alternative Berechnung mit (2*n+1) * Lamda/4 ist auskommentiert im Original)
            if aa1 > _Lamda_viertel:
                val = np.sqrt(aa1 * aa1 - _Lamda_viertel * _Lamda_viertel)
                _Grenzlinie.append(val)
                ij += 1
                Zahl_der_Extrema = ij


# Erzeuge Startpositionen für Feldlinien (x-Koordinaten auf der x-Achse)
_xLinie = []


def LinienStarts_bestimmen():
    global _xLinie
    _xLinie = [0]
    j = 1
    for n in range(Zahl_der_Extrema - 1):
        Linienpos = _Grenzlinie[n] - _Linienabst / 2.0
        iii = 0
        while (Linienpos > 3.0) and (iii < _Linien_pro_Welle):
            _xLinie.append(Linienpos)
            j += 1
            iii += 1
            Linienpos -= _Linienabst
    _xLinie[0] = j - 1
    _xLinie.append(666666.0)  # Endmarkierung (wird nicht genutzt)


def H_Linien_Radien(t):
    """Berechnet alle Radien, bei denen zum Zeitpunkt t magnetische Feldlinien im
    xy-Bild dargestellt werden."""
    radii = []
    r0 = _c * (t - _T / 4.0) - _Wellenlaenge
    max_range = _Wellenlaenge * _Groessenfaktor
    while r0 < max_range:
        for i in range(1, int(_radd[0]) + 1):
            rStart = r0 + _radd[i]
            if rStart > 0 and rStart < max_range:
                radii.append(rStart)
        r0 += _Wellenlaenge / 2.0
    return radii


def E_Linie_berechnen(t, x_start):
    pts = []
    xx = x_start
    yy = 0.0
    Ex0, Ey0, E0 = E_berechnen(t, xx, yy)
    Rechenrichtung = 1 if Ey0 > 0 else -1
    Schritt = 0
    while True:
        rr2 = xx * xx + yy * yy
        if not (
            rr2 >= 9.0 and yy >= 0.0 and Schritt < _ZahlDerRechenschritte and xx > 1.0
        ):
            break
        pts.append((xx, yy))
        xx, yy = Folgepkt(xx, yy, Rechenrichtung, t)
        Schritt += 1
    pts.append((xx, yy))
    Orientierung = -1 if Rechenrichtung > 0 else 1
    return pts, Orientierung


# Runge-Kutta Folgepunkt-Berechnung entlang einer Feldlinie
def Folgepkt(x, y, Rechenrichtung, t):
    Ex0, Ey0, E0 = E_berechnen(t, x, y)
    if E0 == 0:
        return (x, y)
    d = _Grobfaktor2 / E0 * Rechenrichtung
    dx1 = Ex0 * d
    dy1 = Ey0 * d
    Ex1, Ey1, E1 = E_berechnen(t, x + dx1 / 2.0, y + dy1 / 2.0)
    dx2 = Ex1 * d
    dy2 = Ey1 * d
    Ex2, Ey2, E2 = E_berechnen(t, x + dx2 / 2.0, y + dy2 / 2.0)
    dx3 = Ex2 * d
    dy3 = Ey2 * d
    Ex3, Ey3, E3 = E_berechnen(t, x + dx3, y + dy3)
    dx4 = Ex3 * d
    dy4 = Ey3 * d
    x_new = x + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6.0
    y_new = y + (dy1 + 2 * dy2 + 2 * dy3 + dy4) / 6.0
    return (x_new, y_new)


def E_Linie_berechnen2(t, x_start):
    """Integriert eine elektrische Feldlinie (im x-z-Querschnitt) ab Startpunkt (x_start, 0)."""
    pts = []
    xx = x_start
    yy = 0.0
    Ex0, Ey0, _ = E_berechnen(t, xx, yy)
    Rechenrichtung = (
        1 if Ey0 > 0 else -1
    )  # Richtung der Integration (aufwärts oder abwärts)
    Schritt = 0
    while True:
        rr2 = xx * xx + yy * yy
        if not (
            rr2 >= 9.0 and yy >= 0.0 and Schritt < _ZahlDerRechenschritte and xx > 1.0
        ):
            break
        pts.append((xx, yy))
        # Euler-Schritt: gehe einen kleinen Schritt entlang des Feldvektors
        Ex, Ey, _ = E_berechnen(t, xx, yy)
        # Normiere Feldvektor und gehe in dessen Richtung weiter
        E_mag = np.sqrt(Ex * Ex + Ey * Ey)
        if E_mag == 0:
            break
        dx = _Grobfaktor2 * Ex / E_mag
        dy = (_Grobfaktor2 * Ey / E_mag) * Rechenrichtung
        xx += dx
        yy += dy
        Schritt += 1
    pts.append((xx, yy))
    Orientierung = (
        -1 if Rechenrichtung > 0 else 1
    )  # Orientierungsmarkierung (Linie von - zu + oder + zu -)
    return pts, Orientierung


def progress_bar(current, total, length=40):
    percent = int(100 * current / total)
    filled = int(length * current / total)
    bar = "█" * filled + "-" * (length - filled)
    sys.stdout.write(f"\r|{bar}| {percent}% ({current}/{total})")
    sys.stdout.flush()


for frame in range(int(_Periode)):
    # Initiale Berechnungen für t = 0
    _t = (frame) * _dt  # 0.0  # aktuelle Zeit
    progress_bar(frame + 1, _Periode)
    Grenzlinien_bestimmen(_t)
    LinienStarts_bestimmen()

    # Feldlinien (nur obere Halbebene berechnen, Rest durch Spiegelung)
    field_lines_data = []
    field_lines_orient = []
    for i in range(1, int(_xLinie[0]) + 1):
        pts, orient = E_Linie_berechnen(_t, _xLinie[i])
        field_lines_data.append(np.array(pts))
        field_lines_orient.append(orient)

    # Pfeil-Gitter für Energiestrom (nur 1. Quadrant berechnen, Rest später spiegeln)
    xOff = 10.0  # Startversatz auf x-Achse
    x_spacing = _xAbst / _Horizontalfaktor
    y_spacing = _yAbst / _Horizontalfaktor
    # Achsenbereich:
    x_max = _Wellenlaenge * _Groessenfaktor
    y_max = _Wellenlaenge * _Groessenfaktor
    xGrenze = int(round(x_max / (x_spacing)))
    yGrenze = int(round(y_max / (y_spacing)))
    arrow_quadrant_pts = []  # (X, Y) Punkte im 1. Quadranten
    arrow_quadrant_vecs = []  # (Sx, Sy) Werte im 1. Quadranten
    arrow_quadrant_vals = []  # |S| Werte im 1. Quadranten
    for jj in range(yGrenze):
        for ii in range(xGrenze):
            X = xOff + ii * x_spacing
            Y = jj * y_spacing
            Sx, Sy, S = S_berechnen(_t, X, Y)
            arrow_quadrant_pts.append((X, Y))
            arrow_quadrant_vecs.append((Sx, Sy))
            arrow_quadrant_vals.append(S)

    # Aufbau der PyVista-Szene
    plotter = pv.Plotter(off_screen=True)
    plotter.background_color = "white"
    plotter.enable_lightkit(True)  # Beleuchtung für bessere 3D-Wahrnehmung
    plotter.hide_axes()

    # Kamera-Einstellungen für optimale 3D-Ansicht
    radius = (
        _Groessenfaktor * _Wellenlaenge * 1.5
    )  # Größerer Kameraabstand für besseren Überblick
    elev = math.radians(35.0)  # Steilerer Blickwinkel
    azim = math.radians(60.0)  # Anderer Azimut für bessere Sicht auf beide Ebenen
    cam_x = radius * math.cos(elev) * math.cos(azim)
    cam_y = radius * math.cos(elev) * math.sin(azim)
    cam_z = radius * math.sin(elev)
    plotter.camera_position = [(cam_x, cam_y, cam_z), (0, 0, 0), (0, 0, 1)]

    # Anti-Aliasing und Rendering-Qualität verbessern
    plotter.render_window.SetMultiSamples(8)  # Anti-Aliasing
    plotter.renderer.SetUseDepthPeeling(True)  # Bessere Transparenz-Darstellung
    plotter.renderer.SetMaximumNumberOfPeels(8)
    plotter.renderer.SetOcclusionRatio(0.0)

    # Kamerasteuerung und Interaktivität verbessern
    plotter.enable_3_lights()  # Bessere 3D-Beleuchtung
    plotter.camera.zoom(0.8)  # Etwas herauszoomen für besseren Überblick
    plotter.show_axes_all()  # Zeige Orientierungsachsen in der Ecke

    # Laufzeitverhalten und Animation optimieren
    plotter.render_window.SetDesiredUpdateRate(30.0)  # 30 FPS anstreben
    plotter.enable_depth_peeling(10)  # Verbesserte Transparenzdarstellung
    plotter.enable_anti_aliasing()  # Kantenglättung aktivieren

    # Dipolstab und Ladungs-Markierungen
    if Hertzdipol:
        dipole_len = _Lamda_viertel / 8.0  # kurze Dipollänge
    else:
        dipole_len = _Lamda_viertel  # längere Stabdipollänge
    # Dipol-Stab (Linie entlang z-Achse)
    dipolarrow = False
    if dipolarrow:
        arrow_length = (_Lamda_viertel) * (
            -np.sin(_w * _t + 1 / 10)
        )  # Maximale Länge des Pfeils
        dipol_line = pv.Arrow(
            start=(0, 0, -arrow_length / 2),
            direction=(0, 0, arrow_length),
            scale=arrow_length,
            tip_length=0.2,
            shaft_radius=0.05,
        )
        plotter.add_mesh(
            dipol_line,
            color="black",
            name="dipol_arrow",
            line_width=15,
        )
    else:
        dipol_line = pv.Line(
            pointa=(0, 0, -dipole_len / 2), pointb=(0, 0, dipole_len / 2)
        )
        plotter.add_mesh(dipol_line, color="black", line_width=5, name="dipol_line")

    # Plus/Minus-Kugeln an den Dipolenden
    plus_color = "red"
    minus_color = "blue"
    # Kugelradius relativ klein wählen
    sphere_radius = dipole_len * 0.6
    sphere_top = pv.Sphere(center=(0, 0, dipole_len / 2), radius=sphere_radius)
    sphere_bottom = pv.Sphere(center=(0, 0, -dipole_len / 2), radius=sphere_radius)
    sphere_top_actor = plotter.add_mesh(sphere_top, color=plus_color, name="plus_mark")
    sphere_bottom_actor = plotter.add_mesh(
        sphere_bottom, color=minus_color, name="minus_mark"
    )

    # Feldlinien (E-Feld) darstellen
    field_line_polys = []
    field_line_actors = []
    if E_Feldlinien_Flag:
        for pts, orient in zip(field_lines_data, field_lines_orient):
            col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
            X = pts[:, 0]
            Z = pts[:, 1]
            zeros = np.zeros_like(X)
            # Hauptfeldlinie (oben rechts)
            poly_main = pv.lines_from_points(np.column_stack((X, zeros, Z)))
            # Links gespiegelt (oben links)
            poly_left = pv.lines_from_points(np.column_stack((-X, zeros, Z)))
            # Unten gespiegelt (unten rechts)
            poly_bottom = pv.lines_from_points(np.column_stack((X, zeros, -Z)))
            # Unten-links gespiegelt
            poly_bottom_left = pv.lines_from_points(np.column_stack((-X, zeros, -Z)))
            # Alle vier Linien hinzufügen
            field_line_polys.extend(
                [poly_main, poly_left, poly_bottom, poly_bottom_left]
            )
            field_line_actors.extend(
                [
                    plotter.add_mesh(poly_main, color=col, line_width=5),
                    plotter.add_mesh(poly_left, color=col, line_width=5),
                    plotter.add_mesh(poly_bottom, color=col, line_width=5),
                    plotter.add_mesh(poly_bottom_left, color=col, line_width=5),
                ]
            )

    # Magnetfeldlinien (H-Feld) als rote Kreise in der xy-Ebene
    H_line_actors = []
    H_current_radii = []
    if H_Feldlinien_Flag:
        radii_list = H_Linien_Radien(_t)
        for r in radii_list:
            theta = np.linspace(0, 2 * np.pi, 73)
            Xc = r * np.cos(theta)
            Yc = r * np.sin(theta)
            Zc = np.zeros_like(Xc)
            circle = pv.lines_from_points(np.column_stack((Xc, Yc, Zc)), close=True)
            # Magnetfeldlinien in orange-rot mit Transparenz für bessere 3D-Wahrnehmung
            H_line_actors.append(
                plotter.add_mesh(
                    circle,
                    color="#FF4500",  # Orange-Rot opacity=0.7,  # Leichte Transparenz
                    line_width=5,  # Dickere Linien
                    render_lines_as_tubes=True,  # 3D-Röhren statt flacher Linien
                )
            )
        H_current_radii = radii_list.copy()

    # Energiestrom-Pfeile (Poynting-Vektor) via Glyphs
    pfeil_urspruenge = []  # Liste 3D-Punkte aller Pfeilursprünge (alle Quadranten)
    pfeil_richtungen = []  # zugehörige 3D-Richtungsvektoren (werden skaliert berechnet)
    pfeil_magnitude = []  # Betrag des Poynting-Vektors an den Punkten
    # Spiegele Punkte und initiale Vektoren in alle vier Quadranten
    for (X, Y), (Sx, Sy), S in zip(
        arrow_quadrant_pts, arrow_quadrant_vecs, arrow_quadrant_vals
    ):
        # Punkt rechts oben
        pfeil_urspruenge.append([X, Y, 0.0])
        pfeil_richtungen.append([0.0, 0.0, 0.0])  # Platzhalter, wird später gesetzt
        pfeil_magnitude.append(S)
        # Punkt links oben
        pfeil_urspruenge.append([-X, Y, 0.0])
        pfeil_richtungen.append([0.0, 0.0, 0.0])
        pfeil_magnitude.append(S)
        if Y != 0:
            # Punkt rechts unten
            pfeil_urspruenge.append([X, -Y, 0.0])
            pfeil_richtungen.append([0.0, 0.0, 0.0])
            pfeil_magnitude.append(S)
            # Punkt links unten
            pfeil_urspruenge.append([-X, -Y, 0.0])
            pfeil_richtungen.append([0.0, 0.0, 0.0])
            pfeil_magnitude.append(S)

    pfeil_urspruenge = np.array(pfeil_urspruenge)
    pfeil_richtungen = np.array(pfeil_richtungen)
    pfeil_magnitude = np.array(pfeil_magnitude)

    # PolyData für Pfeil-Glyphs erstellen
    pfeil_polydata = pv.PolyData(pfeil_urspruenge)
    # Wir definieren zwei Datenarrays 'vectors' (Richtungsvektor) und 'mag' (Skalierung) auf den Punkten
    pfeil_polydata["vectors"] = pfeil_richtungen
    pfeil_polydata["mag"] = np.zeros_like(
        pfeil_magnitude
    )  # initial noch null, wird in Animation gesetzt
    # Erzeuge einmal einen Pfeil-Quell-Graph (Geometrie) mit gewünschten Proportionen
    arrow_source = pv.Arrow(
        start=(0, 0, 0),
        direction=(1, 0, 0),
        tip_length=0.5,
        tip_radius=0.2,
        shaft_radius=0.04,
    )
    # Initiales Glyph-Objekt (wird aber direkt in Animation erstes Frame aktualisiert)
    pfeil_actor = None
    if Energiestrom_Flag:
        # Glyphs generieren für Pfeile
        initial_glyphs = pfeil_polydata.glyph(
            orient="vectors", scale="mag", factor=1.0, geom=arrow_source
        )
        # Farben ggf. initial einfach alle grau (werden im Update gleich gesetzt)
        pfeil_actor = plotter.add_mesh(initial_glyphs, color="#999999")

    # plotter.reset_camera()
    # Screenshot des Frames erzeugen und zum GIF hinzufügen
    plotter.render()
    plotter.reset_camera_clipping_range()
    # plotter.save_graphic(f"DIPOL{frame}.pdf")
    plotter.screenshot(
        filename=f"./DipolAnimation/DIPOL{frame}.png",
    )
    plotter.close()


# Natural sort key: extract numeric part from filename
def numeric_key(filename):
    match = re.search(r"(\d+)", filename)
    return int(match.group(1)) if match else -1


images = []
file_list = sorted(
    [f for f in os.listdir("DipolAnimation") if f.endswith(".png")], key=numeric_key
)
for filename in file_list:
    images.append(imageio.imread(os.path.join("DipolAnimation", filename)))

imageio.mimsave(
    "EH-field-3D.gif", images, fps=int(round(_Periode / _animation_duration)), loop=0
)
