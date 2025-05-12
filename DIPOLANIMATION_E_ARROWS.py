import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap, BoundaryNorm
import sys

# Physikalische Konstanten und Parameter
_c = 299792458  # Lichtgeschwindigkeit
_c2 = 8.9876e16  # Lichtgeschwindigkeit^2
_epsilon0 = 8.8542e-12
_my0 = 1.2566e-6
_Z0 = 376.7303134

# Dipol-Parameter
_p0 = 100.0  # Amplitude des Dipolvektors
_Wellenlaenge = 256.0  # Wellenlaenge (Simulationseinheiten)
_Lamda_viertel = _Wellenlaenge / 4.0
_w = 2 * np.pi * _c / _Wellenlaenge  # Kreisfrequenz
_Periode = 100  # Anzahl Zeitschritte pro Periode 256 recommended
_animation_duration = 3.0  # Animationsdauer in Sekunden für eine Periode
_T = _Wellenlaenge / _c  # Periodendauer
_dt = _T / _Periode  # Zeitschritt

# Darstellungs- und Steuerungs-Flags
Hertzdipol = True
Stabdipol = False
E_Feldlinien_Flag = True
Energiestrom_Flag = False  # Feldlinien und Energiestrom-Pfeile anzeigen
EnergieMakro_Flag = True  # Makro-Energiestromdichte (wenige, große Pfeile)

# Skalierung Faktoren
_Groessenfaktor = 2
_Horizontalfaktor = 2.0 / _Groessenfaktor
_Vertikalfaktor = 2.0 / _Groessenfaktor

# Parameter für Feldlinien-Berechnung
_Punktzahl_max = 3500
_Grobfaktor2 = 0.25  # Schrittweitenfaktor für Feldlinienintegration
_ZahlDerRechenschritte = int(_Punktzahl_max / _Grobfaktor2 + 100)
_Linienabst = 8.0
_Linien_pro_Welle = _Wellenlaenge / 4.0 / _Linienabst

# Parameter für Energiestrom-Pfeile (Abstände und Pfeilgröße)
if EnergieMakro_Flag:
    _xAbst = 40.0  # Pfeil-Abstand in Pixel (Original)
    _yAbst = 40.0
    _PfeilPktzahl = 20  # Pfeilgrößenfaktor (Original: 20 für Makro)
else:
    _xAbst = 20.0
    _yAbst = 20.0
    _PfeilPktzahl = 10  # Pfeilgrößenfaktor (Original: 10 für Mikro)

# Farbgrenzen für Poynting-Vektor (Energiestromdichte)
mmm = 5.0
Farbgrenze0 = mmm * 1200 * 1e-10
Farbgrenze1 = mmm * 600 * 1e-10
Farbgrenze2 = mmm * 300 * 1e-10
Farbgrenze3 = mmm * 150 * 1e-10
Farbgrenze4 = mmm * 75 * 1e-10
Farbgrenze5 = mmm * 10 * 1e-10

# Farbdefinitionen (aus Original)
Farbe_S0 = "#FF0000"
Farbe_S1 = "#CC0000"
Farbe_S2 = "#AA1111"
Farbe_S3 = "#881188"
Farbe_S4 = "#555555"
Farbe_S5 = "#999999"
Farbe_fieldline_pos = "#008800"  # grün (Linie von -Pol zu +Pol, Orientierung > 0)
Farbe_fieldline_neg = "#0000AA"  # blau (Linie von +Pol zu -Pol, Orientierung < 0)


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


# H-Feld berechnen
def H_berechnen(t, x, y):
    H = 0.0
    if Hertzdipol:
        r2 = x * x + y * y
        r = np.sqrt(r2)
        if r == 0:
            return 0.0
        tt = t - r / _c
        p = _p0 * np.cos(_w * tt)
        pd1 = _w * _p0 * np.sin(_w * tt)  # (Vorzeichen angepasst für z+-Richtung)
        pd2 = _w * _w * p
        sina = x / r
        H = (pd2 / (_c2 * r) * sina) + (pd1 / (_c * r2) * sina)
        return H
    elif Stabdipol:
        rplus = np.sqrt(x * x + (y + _Lamda_viertel) ** 2)
        rminus = np.sqrt(x * x + (y - _Lamda_viertel) ** 2)
        if rplus == 0 or rminus == 0 or x == 0:
            return 0.0
        tt = t - rplus / _c
        sinplus = np.cos(_w * tt)
        tt = t - rminus / _c
        sinminus = np.cos(_w * tt)
        H = (sinplus + sinminus) / (x * 40.0)  # / _Stabfaktor
        return H


# Poynting-Vektor (Energiestromdichte) berechnen
def S_berechnen(t, x, y):
    Ex, Ey, E = E_berechnen(t, x, y)
    H = H_berechnen(t, x, y)
    Sx = Ey * H
    Sy = Ex * H
    S = np.sqrt(Sx * Sx + Sy * Sy)
    return (Sx, Sy, S)


# Energiedichte berechnen (elektrisch + magnetisch)
def Energiedichte_berechnen(t, x, y):
    Ex, Ey, E = E_berechnen(t, x, y)
    H = H_berechnen(t, x, y)
    return _epsilon0 * E * E + _my0 * H * H


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


# Berechnet eine Feldlinie (Liste von (x,y)-Punkten) ab Startpunkt (x_start, 0)
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


# Einrichtung der Grafik
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("equal")
ax.set_xlim(-_Wellenlaenge * _Groessenfaktor, _Wellenlaenge * _Groessenfaktor)
ax.set_ylim(-_Wellenlaenge * _Groessenfaktor, _Wellenlaenge * _Groessenfaktor)
# ax.set_xlabel("x in m")
# ax.set_ylabel("y in m")
ax.set_xlabel(r"x / $\lambda$")
ax.set_ylabel(r"z / $\lambda$")
ax.set_title("Hertzian Dipole - E-Fieldlines")

# Anpassung der Ticks für die x- und y-Achse
ticks = np.arange(
    -_Wellenlaenge * _Groessenfaktor,
    _Wellenlaenge * _Groessenfaktor + 1,
    _Wellenlaenge / 2,
)
ax.set_xticks(ticks)
ax.set_xticklabels(
    [f"{tick / _Wellenlaenge:.2f}" for tick in ticks]
)  # Umrechnung in Meter
ax.set_yticks(ticks)
ax.set_yticklabels(
    [f"{tick / _Wellenlaenge:.2f}" for tick in ticks]
)  # Optional auch für y-Achse

# Initiale Berechnung für t=0
_t = 0.0
Grenzlinien_bestimmen(_t)
LinienStarts_bestimmen()

# Feldlinien-Daten (nur obere Halbebene berechnen, Rest wird durch Spiegelung gezeichnet)
field_lines_data = []
field_lines_orient = []
for i in range(1, int(_xLinie[0]) + 1):
    pts, orient = E_Linie_berechnen(_t, _xLinie[i])
    field_lines_data.append(np.array(pts))
    field_lines_orient.append(orient)

# Pfeil-Gitter (nur 1. Quadrant berechnen, Rest durch Spiegelung zeichnen)
xOff = 10.0  # Startversatz in Simulationseinheiten
x_spacing = _xAbst / _Horizontalfaktor
y_spacing = _yAbst / _Horizontalfaktor
# Bestimme Gittergröße (Anzahl Pfeile) aus Achsenbereich
xGrenze = int(round(ax.get_xlim()[1] / (_xAbst / _Horizontalfaktor)))
yGrenze = int(round(ax.get_ylim()[1] / (_yAbst / _Horizontalfaktor)))
# Energiestrom-Daten (1. Quadrant)
arrow_points = []
arrow_vectors = []
arrow_vals = []
for jj in range(yGrenze):
    for ii in range(xGrenze):
        X = xOff + ii * x_spacing
        Y = jj * y_spacing
        Sx, Sy, S = S_berechnen(_t, X, Y)
        arrow_points.append((X, Y))
        arrow_vectors.append((Sx, Sy))
        arrow_vals.append(S)
# Spiegelung der Pfeile in alle vier Quadranten
full_X = []
full_Y = []
full_U = []
full_V = []
full_S = []
for (X, Y), (Sx, Sy), S in zip(arrow_points, arrow_vectors, arrow_vals):
    # rechts oben
    full_X.append(X)
    full_Y.append(Y)
    full_U.append(0.0)
    full_V.append(0.0)
    full_S.append(S)  # Platzhalter, wird unten berechnet
    # links oben (X negiert)
    full_X.append(-X)
    full_Y.append(Y)
    full_U.append(0.0)
    full_V.append(0.0)
    full_S.append(S)
    if Y != 0:
        # rechts unten (Y negiert)
        full_X.append(X)
        full_Y.append(-Y)
        full_U.append(0.0)
        full_V.append(0.0)
        full_S.append(S)
        # links unten (X und Y negiert)
        full_X.append(-X)
        full_Y.append(-Y)
        full_U.append(0.0)
        full_V.append(0.0)
        full_S.append(S)
# Normalisiere Pfeil-Vektoren auf fixe Länge (Skalierung analog Original: 1.2 * _PfeilPktzahl)
full_X = np.array(full_X)
full_Y = np.array(full_Y)
full_Sx = []
full_Sy = []
for (X, Y), (Sx, Sy) in zip(arrow_points, arrow_vectors):
    # Vektoren aus 1. Quadrant normalisieren und später spiegeln
    if Sx == 0 and Sy == 0:
        full_Sx.append(0.0)
        full_Sy.append(0.0)
    else:
        S_val = np.sqrt(Sx * Sx + Sy * Sy)
        scale = 1.2 * _PfeilPktzahl / (S_val if S_val != 0 else 1.0)
        full_Sx.append(Sx * scale)
        full_Sy.append(Sy * scale)
# Jetzt Spiegelwerte zuordnen:
# full_U, full_V in gleicher Reihenfolge wie in for-Schleife initialisiert
k = 0
for (X, Y), (Sx_scaled, Sy_scaled) in zip(arrow_points, zip(full_Sx, full_Sy)):
    U_val, V_val = Sx_scaled, Sy_scaled
    # rechts oben
    full_U[k] = U_val
    full_V[k] = V_val
    k += 1
    # links oben
    full_U[k] = -U_val
    full_V[k] = V_val
    k += 1
    if Y != 0:
        # rechts unten
        full_U[k] = U_val
        full_V[k] = -V_val
        k += 1
        # links unten
        full_U[k] = -U_val
        full_V[k] = -V_val
        k += 1

full_U = np.array(full_U)
full_V = np.array(full_V)
full_S = np.array(full_S)

# Colormap für Energiestrom-Pfeile (diskrete Farben)
colors_list = [Farbe_S5, Farbe_S4, Farbe_S3, Farbe_S2, Farbe_S1, Farbe_S0]
boundaries = [
    0.0,
    Farbgrenze5,
    Farbgrenze4,
    Farbgrenze3,
    Farbgrenze2,
    Farbgrenze1,
    np.inf,
]
cmap = ListedColormap(colors_list)
norm = BoundaryNorm(boundaries, cmap.N)

# Zeichne anfängliche Feldlinien (und gespiegelte)
field_line_plots = []  # Liste für Feldlinien-Objekte
arrow_objs = []  # Liste für Pfeil-Objekte
for pts, orient in zip(field_lines_data, field_lines_orient):
    col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
    X = pts[:, 0]
    Y = pts[:, 1]
    # Hauptfeldlinie (obere rechte)
    (ln_main,) = ax.plot(X, Y, color=col, linewidth=1)
    # Horizontal gespiegelte Feldlinie (obere linke)
    (ln_left,) = ax.plot(-X, Y, color=col, linewidth=1)
    # Vertikal gespiegelte Feldlinie (untere rechte)
    (ln_bottom,) = ax.plot(X, -Y, color=col, linewidth=1)
    # Doppelt gespiegelte Feldlinie (untere linke)
    (ln_bottom_left,) = ax.plot(-X, -Y, color=col, linewidth=1)
    field_line_plots.extend([ln_main, ln_left, ln_bottom, ln_bottom_left])

if E_Feldlinien_Flag:  # nur wenn Linien überhaupt gezeigt werden
    for pts, orient in zip(field_lines_data, field_lines_orient):
        col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
        mark = "^" if orient > 0 else "v"  # ^ = Pfeil nach oben, v = nach unten
        x0 = pts[0, 0]  # Startpunkt dieser Feldlinie (y = 0)
        # rechts
        (arr_r,) = ax.plot(
            x0, 0, marker=mark, markersize=4, color=col, linestyle="None"
        )
        # links (gespiegelte Linie)
        (arr_l,) = ax.plot(
            -x0, 0, marker=mark, markersize=4, color=col, linestyle="None"
        )
        arrow_objs.extend([arr_r, arr_l])  # merken fürs spätere Update

# Zeichne anfängliche Energiestrom-Pfeile
# quiver = ax.quiver(
#     full_X,
#     full_Y,
#     full_U,
#     full_V,
#     full_S,
#     cmap=cmap,
#     norm=norm,
#     angles="xy",
#     scale_units="xy",
#     scale=1,
#     width=0.004,
# )

# Zeichne Dipol (schwarzer Strich) und Markierung für +/-
dipole_objs = []
if Hertzdipol:
    dipole_len = _Lamda_viertel / 8.0  # kleiner Strich
    (dipole_line,) = ax.plot(
        [0, 0], [-dipole_len / 2, dipole_len / 2], color="black", linewidth=2, zorder=3
    )
    dipole_objs.append(dipole_line)
elif Stabdipol:
    dipole_len = _Lamda_viertel  # halbe Dipollänge (Lamda/4 in jede Richtung)
    (dipole_line,) = ax.plot(
        [0, 0], [-dipole_len, dipole_len], color="black", linewidth=2, zorder=3
    )
    dipole_objs.append(dipole_line)
# Plus/Minus-Markierung (rot = +, blau = -)
p_val = _p0 * np.cos(_w * _t)
if Hertzdipol:
    if p_val > 0:
        # oben + (rot), unten - (blau)
        (plus_mark,) = ax.plot(0, dipole_len / 2, "o", color="red", zorder=4)
        (minus_mark,) = ax.plot(0, -dipole_len / 2, "o", color="blue", zorder=4)
    else:
        (plus_mark,) = ax.plot(0, -dipole_len / 2, "o", color="red", zorder=4)
        (minus_mark,) = ax.plot(0, dipole_len / 2, "o", color="blue", zorder=4)
    dipole_objs.extend([plus_mark, minus_mark])
elif Stabdipol:
    if p_val > 0:
        (plus_mark,) = ax.plot(0, dipole_len, "o", color="red", zorder=4)
        (minus_mark,) = ax.plot(0, -dipole_len, "o", color="blue", zorder=4)
    else:
        (plus_mark,) = ax.plot(0, -dipole_len, "o", color="red", zorder=4)
        (minus_mark,) = ax.plot(0, dipole_len, "o", color="blue", zorder=4)
    dipole_objs.extend([plus_mark, minus_mark])


def progress_bar(current, total, length=40):
    percent = int(100 * current / total)
    filled = int(length * current / total)
    bar = "█" * filled + "-" * (length - filled)
    sys.stdout.write(f"\r|{bar}| {percent}% ({current}/{total})")
    sys.stdout.flush()


# Animations-Update-Funktion
def update(frame):
    # print(f"Frame: {frame + 1} of {_Periode}")
    progress_bar(frame + 1, _Periode)
    global _t, arrow_objs
    _t = frame * _dt
    if _t > _T:
        _t -= _T  # Zeitsprung (Periodisch)
    # Feldlinien neu berechnen
    Grenzlinien_bestimmen(_t)
    LinienStarts_bestimmen()
    new_field_lines = []
    new_orientations = []
    for i in range(1, int(_xLinie[0]) + 1):
        pts, orient = E_Linie_berechnen(_t, _xLinie[i])
        new_field_lines.append(np.array(pts))
        new_orientations.append(orient)
    # Passen wir ggf. die Anzahl der gezeichneten Linien an
    total_segments = len(new_field_lines) * 4
    if total_segments != len(field_line_plots):
        # Entferne alte und zeichne neue (Anzahl hat sich geändert)
        for ln in field_line_plots:
            ln.remove()
        field_line_plots.clear()
        for pts, orient in zip(new_field_lines, new_orientations):
            col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
            X = pts[:, 0]
            Y = pts[:, 1]
            (ln_main,) = ax.plot(X, Y, color=col, linewidth=1)
            (ln_left,) = ax.plot(-X, Y, color=col, linewidth=1)
            (ln_bottom,) = ax.plot(X, -Y, color=col, linewidth=1)
            (ln_bottom_left,) = ax.plot(-X, -Y, color=col, linewidth=1)
            field_line_plots.extend([ln_main, ln_left, ln_bottom, ln_bottom_left])
        # ---------- alte Pfeile löschen, neue anlegen ----------
        for ar in arrow_objs:
            ar.remove()
        arrow_objs = []
        if E_Feldlinien_Flag:
            for pts, orient in zip(new_field_lines, new_orientations):
                col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
                mark = "^" if orient > 0 else "v"
                x0 = pts[0][0]
                (arr_r,) = ax.plot(
                    x0, 0, marker=mark, color=col, markersize=4, linestyle="None"
                )
                (arr_l,) = ax.plot(
                    -x0, 0, marker=mark, color=col, markersize=4, linestyle="None"
                )
                arrow_objs.extend([arr_r, arr_l])
        # --------------------------------------------------------
    else:
        # Aktualisiere bestehende Linien
        idx = 0
        for pts, orient in zip(new_field_lines, new_orientations):
            col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
            X = pts[:, 0]
            Y = pts[:, 1]
            field_line_plots[idx].set_data(X, Y)
            field_line_plots[idx].set_color(col)
            idx += 1
            field_line_plots[idx].set_data(-X, Y)
            field_line_plots[idx].set_color(col)
            idx += 1
            field_line_plots[idx].set_data(X, -Y)
            field_line_plots[idx].set_color(col)
            idx += 1
            field_line_plots[idx].set_data(-X, -Y)
            field_line_plots[idx].set_color(col)
            idx += 1
        # ---------- Pfeilspitzen nur versetzen / einfärben ----------
        if E_Feldlinien_Flag:
            idx_arr = 0
            for pts, orient in zip(new_field_lines, new_orientations):
                col = Farbe_fieldline_pos if orient > 0 else Farbe_fieldline_neg
                mark = "^" if orient > 0 else "v"
                x0 = pts[0][0]

                arrow_objs[idx_arr].set_data([x0], [0])
                arrow_objs[idx_arr].set_marker(mark)
                arrow_objs[idx_arr].set_color(col)
                idx_arr += 1

                arrow_objs[idx_arr].set_data([-x0], [0])
                arrow_objs[idx_arr].set_marker(mark)
                arrow_objs[idx_arr].set_color(col)
                idx_arr += 1
        else:
            for ar in arrow_objs:
                ar.set_data([], [])
        # ------------------------------------------------------------
    # Berechne Energiestrom-Pfeile neu
    new_U = []
    new_V = []
    new_Svals = []
    # Erste Quadrant berechnen und gleich normalisieren
    for (X, Y), (Sx, Sy), S in zip(arrow_points, arrow_vectors, arrow_vals):
        # Neu berechnen an ursprünglichen Gitterpunkten
        Sx_new, Sy_new, S_new = S_berechnen(_t, X, Y)
        # Normalisieren
        if S_new == 0:
            U_val = 0.0
            V_val = 0.0
        else:
            scale = 1.2 * _PfeilPktzahl / S_new
            U_val = Sx_new * scale
            V_val = Sy_new * scale
        # Spiegeln in alle Quadranten
        new_U.append(U_val)
        new_V.append(V_val)
        new_Svals.append(S_new)
        new_U.append(-U_val)
        new_V.append(V_val)
        new_Svals.append(S_new)
        if Y != 0:
            new_U.append(U_val)
            new_V.append(-V_val)
            new_Svals.append(S_new)
            new_U.append(-U_val)
            new_V.append(-V_val)
            new_Svals.append(S_new)
    new_U = np.array(new_U)
    new_V = np.array(new_V)
    new_Svals = np.array(new_Svals)
    # Update Quiver (Richtungen und Farben)
    # quiver.set_UVC(new_U, new_V, new_Svals)

    # Aktualisiere Dipol-Markierungen (+/-)
    p_val = _p0 * np.cos(_w * _t)
    if Hertzdipol:
        if p_val > 0:
            dipole_objs[-2].set_data([0], [dipole_len / 2])
            dipole_objs[-2].set_color("blue")
            dipole_objs[-1].set_data([0], [-dipole_len / 2])
            dipole_objs[-1].set_color("red")
        else:
            dipole_objs[-2].set_data([0], [-dipole_len / 2])
            dipole_objs[-2].set_color("blue")
            dipole_objs[-1].set_data([0], [dipole_len / 2])
            dipole_objs[-1].set_color("red")
    elif Stabdipol:
        if p_val > 0:
            dipole_objs[-2].set_data([0], [dipole_len])
            dipole_objs[-2].set_color("blue")
            dipole_objs[-1].set_data([0], [-dipole_len])
            dipole_objs[-1].set_color("red")
        else:
            dipole_objs[-2].set_data([0], [-dipole_len])
            dipole_objs[-2].set_color("blue")
            dipole_objs[-1].set_data([0], [dipole_len])
            dipole_objs[-1].set_color("red")
    return tuple(field_line_plots) + tuple(dipole_objs) + tuple(arrow_objs)


# Erstelle Animation
anim = FuncAnimation(fig, update, frames=int(_Periode), interval=100, blit=True)
# Animation automatisch als GIF-Datei (eine Periodendauer) speichern
print("\n Creating gif...")
anim.save(
    "E-field-2D.gif",
    writer="pillow",
    fps=int(round(_Periode / _animation_duration)),
    dpi=300,
)
print("\nAnimation saved as 'dipol_E_field.gif'")

# plt.show()
