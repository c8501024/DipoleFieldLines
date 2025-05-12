import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Bbox
from matplotlib.patches import Circle


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
_Periode = 100  # Anzahl Zeitschritte pro Periode
_T = _Wellenlaenge / _c  # Periodendauer
_dt = _T / _Periode  # Zeitschritt
_animation_duration = 3  # Dauer der Animation in Sekunden

# Darstellungs- und Steuerungs-Flags
Hertzdipol = True
Stabdipol = False
E_Feldlinien_Flag = True
H_Feldlinien_Flag = True  # magnetische Feldlinien (H-Feld) anzeigen
Energiestrom_Flag = False  # Feldlinien und Energiestrom-Pfeile anzeigen
EnergieMakro_Flag = True  # Makro-Energiestromdichte (wenige, große Pfeile)

# Skalierung Faktoren
_Groessenfaktor = 1.5
_Horizontalfaktor = 2.0 / _Groessenfaktor
_Vertikalfaktor = 2.0 / _Groessenfaktor

# Parameter für Feldlinien-Berechnung
_Punktzahl_max = 3500
_Grobfaktor2 = 0.25  # Schrittweitenfaktor für Feldlinienintegration
_ZahlDerRechenschritte = int(_Punktzahl_max / _Grobfaktor2 + 100)
_Linienabst = 8.0
_Linien_pro_Welle = _Wellenlaenge / 4.0 / _Linienabst

# Parameter für magnetische Feldlinien in der xy-Ebene (H-Feld)
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


# **Hilfsfunktion zur Berechnung der Radien der H-Feldlinien in der xy-Ebene**
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
            x2 = x1 + Schrittweite
            if x1 > _Wellenlaenge * _Groessenfaktor:
                Fensterende = True


# Runge-Kutta Folgepunkt-Berechnung entlang einer Feldlinie
def Folgepkt(x, y, Rechenrichtung, t):
    """
    Berechnet den nächsten Punkt (x_new, y_new) entlang einer Feldlinie
    mittels RK4-Integration.  Schrittweite ~ 1/E  * _Grobfaktor2.
    Rechenrichtung = +1  ⇒ mit dem Feld
                    = -1 ⇒ entgegen dem Feld
    """
    Ex0, Ey0, E0 = E_berechnen(t, x, y)
    if E0 == 0:
        return x, y
    if E0 < 1e-6:  # Sicherheit, sonst riesiger Schritt
        E0 = 1e-6
    h = _Grobfaktor2 / E0 * Rechenrichtung  # Schrittweite

    dx1, dy1 = Ex0 * h, Ey0 * h
    Ex1, Ey1, _ = E_berechnen(t, x + dx1 / 2, y + dy1 / 2)
    dx2, dy2 = Ex1 * h, Ey1 * h
    Ex2, Ey2, _ = E_berechnen(t, x + dx2 / 2, y + dy2 / 2)
    dx3, dy3 = Ex2 * h, Ey2 * h
    Ex3, Ey3, _ = E_berechnen(t, x + dx3, y + dy3)
    dx4, dy4 = Ex3 * h, Ey3 * h

    x_new = x + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6
    y_new = y + (dy1 + 2 * dy2 + 2 * dy3 + dy4) / 6
    return x_new, y_new


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


# Einrichtung der Grafik
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)  # 3D-Achse

ax.set_aspect("equal")
ax.set_xlim(-_Wellenlaenge * _Groessenfaktor, _Wellenlaenge * _Groessenfaktor)
ax.set_ylim(-_Wellenlaenge * _Groessenfaktor, _Wellenlaenge * _Groessenfaktor)


ax.set_xlabel(r"x / $\lambda$")
ax.set_ylabel(r"y / $\lambda$")

ax.set_title("Hertzian Dipole - H-Fieldlines")

# Anpassung der Ticks für die x- und y-Achse
ticks = np.arange(
    -_Wellenlaenge * _Groessenfaktor,
    _Wellenlaenge * _Groessenfaktor + 1,
    _Wellenlaenge / 2,
)
ax.set_xticks(ticks)
ax.set_xticklabels([f"{tick / _Wellenlaenge:.2f}" for tick in ticks])
ax.set_yticks(ticks)
ax.set_yticklabels([f"{tick / _Wellenlaenge:.2f}" for tick in ticks])

# Initiale Berechnung für t=0
_t = 0.0
Grenzlinien_bestimmen(_t)
LinienStarts_bestimmen()

# Feldlinien-Daten (nur obere Halbebene berechnen, Rest wird durch Spiegelung gezeichnet)


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

# Zeichne anfängliche Feldlinien (und gespiegelte)

mag_arrow_objects = []  # Liste für Pfeil-Objekte
# **Zeichne magnetische Feldlinien (rot) in der xy-Ebene**
H_line_plots = []  # Liste für H-Feldlinien-Objekte
H_current_radii = []  # aktuelle Radien der H-Feldlinien
if H_Feldlinien_Flag:
    radii_list = H_Linien_Radien(_t)
    for r in radii_list:
        theta = np.linspace(0, 2 * np.pi, 73)
        Xc = r * np.cos(theta)
        Yc = r * np.sin(theta)
        (ln_h,) = ax.plot(Xc, Yc, color="red", linewidth=1, zorder=2)
        H_line_plots.append(ln_h)
    H_current_radii = radii_list.copy()

    for r in radii_list:
        Hphi = H_berechnen(0, r, 0)
        mark = "^" if Hphi > 0 else "v"
        col = "red" if Hphi > 0 else "blue"

        x0 = r  # Startpunkt dieser Feldlinie (y = 0)
        # rechts
        (arr_r,) = ax.plot(
            x0, 0, marker=mark, markersize=6, color=col, linestyle="None"
        )
        # links (gespiegelte Linie)
        (arr_l,) = ax.plot(
            -x0, 0, marker=mark, markersize=6, color=col, linestyle="None"
        )
        mag_arrow_objects.extend([arr_r, arr_l])  # merken fürs spätere Update

# Zeichne Dipol (schwarzer Strich) und Markierung für +/-
dipole_objs = []
dipole_objs = []

# dipol zeichnen
# Erstelle den Kreis und das Kreuz einmalig
# circle = Circle((0, 0), radius=2, edgecolor="black", facecolor="white", zorder=5)
# ax.add_patch(circle)

# Linien für das Kreuz
(cross_line1,) = ax.plot([], [], "kx", markersize=10, zorder=6)  # Diagonale 1
# Punkt für den Kreis
(circle_point,) = ax.plot([], [], "ko", markersize=10, zorder=6)  # Schwarzer Punkt


def progress_bar(current, total, length=40):
    percent = int(100 * current / total)
    filled = int(length * current / total)
    bar = "█" * filled + "-" * (length - filled)
    sys.stdout.write(f"\r|{bar}| {percent}% ({current}/{total})")
    sys.stdout.flush()


def update(frame):
    global _t, mag_arrow_objects
    _t = frame * _dt
    if _t > _T:
        _t -= _T  # Zeitsprung (periodisch)

    progress_bar(frame + 1, _Periode)
    # Feldlinien neu berechnen
    Grenzlinien_bestimmen(_t)
    LinienStarts_bestimmen()

    # **Magnetfeldlinien (H-Feld) neu berechnen und aktualisieren**
    if H_Feldlinien_Flag:
        new_radii = H_Linien_Radien(_t)
        if (len(new_radii) != len(H_current_radii)) or (
            sorted(new_radii) != sorted(H_current_radii)
        ):
            # Anzahl oder Position der H-Feldlinien hat sich geändert -> neu zeichnen
            for ln in H_line_plots:
                ln.remove()
            H_line_plots.clear()
            for r in new_radii:
                theta = np.linspace(0, 2 * np.pi, 73)
                Xc = r * np.cos(theta)
                Yc = r * np.sin(theta)
                (ln_h,) = ax.plot(Xc, Yc, color="red", linewidth=1, zorder=2)
                H_line_plots.append(ln_h)

            H_current_radii.clear()
            H_current_radii.extend(new_radii)
        else:
            # Positionen der bestehenden H-Feldlinien aktualisieren
            for i, r in enumerate(new_radii):
                theta = np.linspace(0, 2 * np.pi, 73)
                Xc = r * np.cos(theta)
                Yc = r * np.sin(theta)
                H_line_plots[i].set_data(Xc, Yc)
    else:
        # Magnetfeldlinien verbergen (keine Anzeige)
        for ln in H_line_plots:
            ln.remove()
        H_line_plots.clear()
        H_current_radii.clear()

    for ar in mag_arrow_objects:
        ar.remove()
    mag_arrow_objects = []

    for r in new_radii:
        Hphi = H_berechnen(_t, r, 0)
        mark1 = "^" if Hphi > 0 else "v"
        mark2 = "v" if Hphi > 0 else "^"
        col = "red" if Hphi > 0 else "blue"

        x0 = r  # Startpunkt dieser Feldlinie (y = 0)
        # rechts
        (arr_r,) = ax.plot(
            x0, 0, marker=mark1, markersize=6, color=col, linestyle="None"
        )
        # links (gespiegelte Linie)
        (arr_l,) = ax.plot(
            -x0, 0, marker=mark2, markersize=6, color=col, linestyle="None"
        )
        mag_arrow_objects.extend([arr_r, arr_l])  # merken fürs spätere Update
    # Aktualisiere Dipol-Markierungen (+/-) und zeichne Kreis oder Kreuz

    Hphi = H_berechnen(_t, 0.1, 0)

    if Hphi > 0:
        # Kreis mit Punkt (Strom in +z-Richtung)
        circle_point.set_data([0], [0])  # Punkt in der Mitte
        cross_line1.set_data([], [])  # Kreuz unsichtbar
    else:
        # Kreis mit Kreuz (Strom in -z-Richtung)
        circle_point.set_data([], [])  # Punkt unsichtbar
        cross_line1.set_data([-0.35, 0.35], [-0.35, 0.35])  # Diagonale 1

    return tuple(H_line_plots) + (circle_point, cross_line1) + tuple(mag_arrow_objects)


# Erstelle Animation
anim = FuncAnimation(fig, update, frames=int(_Periode), interval=100, blit=False)

# Animation automatisch als GIF-Datei speichern (optional)
print("\n Creating gif...")
anim.save(
    "H-field_2D.gif",
    writer="pillow",
    fps=int(round(_Periode / _animation_duration)),
    dpi=300,
)
print("\nAnimation saved as 'dipol_E_field.gif'")
# plt.show()
