# -*- coding: utf-8 -*-
"""
Animation der E-Feldlinien eines Halbwellendipols (Folie 19 / 21 der Präsentation).

Grundlage ist die exakte (geschlossene) Lösung für die sinusförmige Stromverteilung
I(z) = I0 cos(kz) auf einem dünnen Dipol der Länge L = lambda/2 (siehe Folie 21):

    H_phi = j * I0 * e^{j w t} / (4 pi rho) * ( e^{-j k R1} + e^{-j k R2} )
    R_{1,2} = sqrt( rho^2 + (z -/+ lambda/4)^2 )      (Abstände zu den Dipolenden)

Feldlinien: Bei einem rotationssymmetrischen Feld (nur H_phi, E_rho, E_z) folgt aus rot H = j w eps0 E

    E = 1/(j w eps0 rho) * ( -d(rho H_phi)/dz ,  d(rho H_phi)/drho ),

d. h. E liegt tangential zu den Höhenlinien von rho * H_phi. Für den Momentanwert Re{E e^{j w t}}
wandert der Faktor 1/j = -j mit in die Höhenlinienfunktion; die momentanen E-Feldlinien sind daher
die Höhenlinien von

    psi(rho, z, t) = Re{ -j * rho * H_phi * e^{j w t} } ~ Re{ (e^{-j k R1} + e^{-j k R2}) e^{j w t} }

Der Faktor -j bedeutet eine Viertelperiode Verzögerung gegenüber H (E entsteht durch zeitliche
Integration von rot H): Beim Strommaximum ist das Nahfeld von E schwach, bei maximaler
Ladungstrennung (Strom = 0) am stärksten. Konstante Vorfaktoren spielen für die Form keine Rolle.
Hinweis: Die Liniendichte ist proportional zu rho*|E|, nicht exakt zu |E|.
Dipol: schwarzer Stab der Länge lambda/2 mit zwei stellvertretenden Ladungsträgern +q (rot) und -q (blau).
Sie schwingen gegenläufig:  z_+(t) = +h sin(wt),  z_-(t) = -h sin(wt)   (h = lambda/4).
Ihre Geschwindigkeit ist ~ cos(wt) ~ I(t) = I0 cos(wt): Beim Strommaximum laufen sie mit größter
Geschwindigkeit durch die Mitte, beim Stromnulldurchgang stehen sie an den Dipolenden (maximale
Ladungstrennung, wie die Ladungskugeln beim Hertzschen Dipol). +q nach oben und -q nach unten
entsprechen beide einem Strom in +z-Richtung.

Benötigt: numpy, matplotlib, pillow
"""
import io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from PIL import Image

# ------------------------------------------------------------------ Parameter
lam = 1.0                 # Wellenlänge (alle Längen in Vielfachen von lambda)
k = 2 * np.pi / lam       # Wellenzahl
h = lam / 4               # halbe Dipollänge (L = lambda/2)
N_FRAMES = 80             # Bilder pro Periode
FRAME_MS = 30             # Anzeigedauer je Bild in ms (30 ms ~ 33 Bilder/s)
LW = 1.5                  # Linienstärke der Feldlinien
N_GRID = 601              # Auflösung des Rechengitters
EXTENT = 2.0              # dargestellter Bereich: -2 lambda ... +2 lambda
LEVELS = np.linspace(-1.9, 1.9, 20)   # Höhenlinien (psi liegt zwischen -2 und 2)
FARBE_PLUS, FARBE_MINUS = 'red', 'blue'   # Ladungsträger +q und -q
MS = 11                   # Größe der Ladungsträger (Punkte)
OUT = 'feldlinien_halbwellendipol.gif'

# ------------------------------------------------------------------ Feld auf dem Gitter
x = np.linspace(-EXTENT, EXTENT, N_GRID)
z = np.linspace(-EXTENT, EXTENT, N_GRID)
X, Z = np.meshgrid(x, z)
rho = np.abs(X) + 1e-9                         # Abstand zur Dipolachse
R1 = np.sqrt(rho**2 + (Z - h)**2)              # Abstand zum oberen Dipolende
R2 = np.sqrt(rho**2 + (Z + h)**2)              # Abstand zum unteren Dipolende
P = np.exp(-1j * k * R1) + np.exp(-1j * k * R2)   # -j*rho*H_phi bis auf den Faktor I0/(4 pi)

# ------------------------------------------------------------------ Bilder erzeugen
frames = []
for i in range(N_FRAMES):
    wt = 2 * np.pi * i / N_FRAMES                  # Phase omega*t
    psi = np.real(P * np.exp(1j * wt))             # Höhenlinienfunktion der momentanen E-Feldlinien

    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    ax.contour(X, Z, psi, levels=LEVELS[LEVELS > 0], colors='green', linewidths=LW, linestyles='solid')
    ax.contour(X, Z, psi, levels=LEVELS[LEVELS < 0], colors='darkblue', linewidths=LW, linestyles='solid')

    # Dipol: schwarzer Stab mit zwei gegenläufig schwingenden Ladungsträgern (Geschwindigkeit ~ I(t))
    ax.plot([0, 0], [-h, h], color='k', lw=4, solid_capstyle='butt', zorder=5)
    z_plus = h * np.sin(wt)
    for zq, farbe, zeichen, zo in [(-z_plus, FARBE_MINUS, '\u2212', 6), (z_plus, FARBE_PLUS, '+', 7)]:
        ax.plot(0, zq, 'o', color=farbe, ms=MS, mec='white', mew=0.8, zorder=zo)
        ax.text(0, zq, zeichen, color='white', fontsize=MS * 0.75, fontweight='bold',
                ha='center', va='center_baseline', zorder=zo)

    ax.set_aspect('equal')
    ax.set_xlim(-EXTENT, EXTENT); ax.set_ylim(-EXTENT, EXTENT)
    ax.set_xticks(np.arange(-EXTENT, EXTENT + 0.01, 0.5)); ax.set_yticks(np.arange(-EXTENT, EXTENT + 0.01, 0.5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f')); ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel(r'x / $\lambda$'); ax.set_ylabel(r'z / $\lambda$')
    ax.set_title('Halbwellendipol - Feldlinien')

    buf = io.BytesIO()
    fig.savefig(buf, format='png', facecolor='white'); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert('RGB'))

# ------------------------------------------------------------------ weißen Rand abschneiden, GIF speichern
arrs = [np.asarray(f) for f in frames]
mask = np.any(np.stack([(a < 245).any(axis=2) for a in arrs]), axis=0)
ys, xs = np.where(mask)
box = (max(xs.min() - 20, 0), max(ys.min() - 20, 0),
       min(xs.max() + 20, arrs[0].shape[1]), min(ys.max() + 20, arrs[0].shape[0]))
frames = [f.crop(box) for f in frames]
# gemeinsame Farbpalette für alle Bilder (sonst gehen seltene Farben wie die der Ladungsträger verloren)
w, h = frames[0].size
sheet = Image.new('RGB', (w, h * len(frames)))
for j, f in enumerate(frames):
    sheet.paste(f, (0, j * h))
base = sheet.resize((w // 2, h * len(frames) // 2)).quantize(colors=26, method=Image.Quantize.MEDIANCUT)
fixed = [(255, 255, 255), (0, 0, 0), (255, 0, 0), (0, 0, 255), (0, 128, 0), (0, 0, 139), (62, 48, 103)]   # Farben der Grafik
pal_rgb = np.array(base.getpalette()[:26 * 3] + [c for rgb in fixed for c in rgb], dtype=np.int32).reshape(-1, 3)


def auf_palette(f):
    """Jedes Pixel exakt der nächstgelegenen Palettenfarbe zuordnen (Weiß bleibt reinweiß)."""
    a = np.asarray(f).reshape(-1, 3)
    cols, inv = np.unique(a, axis=0, return_inverse=True)
    d = ((cols[:, None, :].astype(np.int32) - pal_rgb[None, :, :]) ** 2).sum(-1)
    idx = d.argmin(1).astype(np.uint8)[inv.ravel()]
    im = Image.fromarray(idx.reshape(f.size[1], f.size[0]), 'P')
    im.putpalette(pal_rgb.astype(np.uint8).ravel().tolist())
    return im


pal = [auf_palette(f) for f in frames]
pal[0].save(OUT, save_all=True, append_images=pal[1:], duration=FRAME_MS, loop=0, optimize=True)
print('gespeichert:', OUT, frames[0].size)
