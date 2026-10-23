# -*- coding: utf-8 -*-
"""
Animation der E-Feldlinien eines Hertzschen Dipols (Gegenstück zu feldlinien_halbwellendipol.py).

Grundlage ist das exakte Feld des Hertzschen Dipols (Länge d << lambda, konstanter Strom I0),
wie in der Vorlesung hergeleitet:

    H_phi = k^2 I0 d / (4 pi) * e^{j w t} e^{-j k r} * ( j/(kr) + 1/(kr)^2 ) * sin(theta)

Feldlinien: Bei einem rotationssymmetrischen Feld (nur H_phi, E_r, E_theta) folgt aus
rot H = j w eps0 E, dass E überall tangential zu den Höhenlinien von  rho * H_phi  liegt
(rho = r sin(theta) = Abstand zur Dipolachse). Für den Momentanwert Re{E e^{j w t}} wandert der
Faktor 1/(j w eps0) mit in die Höhenlinienfunktion. Mit rho * H_phi ~ sin^2(theta) (j + 1/(kr)) e^{-jkr}
sind die momentanen E-Feldlinien daher die Höhenlinien von

    psi(r, theta, t) = Re{ -j * rho * H_phi * e^{j w t} } ~ Re{ sin^2(theta) (1 - j/(kr)) e^{j(w t - k r)} }
                     = sin^2(theta) * [ cos(w t - k r) + sin(w t - k r)/(k r) ]

Im Nahfeld ist psi ~ sin^2(theta) sin(w t)/(kr) ~ q(t): die quasistatischen Feldlinien eines
Dipols mit dem momentanen Moment p(t) = q(t) d, passend zu den Ladungskugeln.
Konstante Vorfaktoren spielen für die Form der Linien keine Rolle.
Hinweis: Die Liniendichte ist proportional zu rho*|E|, nicht exakt zu |E|.
Dipol: kurzer schwarzer Stab mit zwei stellvertretenden Ladungsträgern +q (rot) und -q (blau), die
gegenläufig schwingen:  z_+(t) = +(d/2) sin(wt),  z_-(t) = -(d/2) sin(wt).
Das Dipolmoment p(t) ~ sin(wt) und der Strom I(t) = dp/dt / d ~ cos(wt) passen damit zur Rechnung:
Beim Strommaximum laufen die Ladungen durch die Mitte (p = 0), beim Stromnulldurchgang sind sie
maximal getrennt. Die Stablänge d_plot dient nur der Darstellung (gerechnet wird mit d << lambda).

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
d_plot = 0.24             # nur für die Darstellung: Stablänge = maximaler Abstand der Ladungen (in lambda)
N_FRAMES = 80             # Bilder pro Periode
FRAME_MS = 30             # Anzeigedauer je Bild in ms (30 ms ~ 33 Bilder/s)
LW = 1.5                  # Linienstärke der Feldlinien
N_GRID = 601              # Auflösung des Rechengitters
EXTENT = 2.0              # dargestellter Bereich: -2 lambda ... +2 lambda
LEVELS = np.linspace(-0.95, 0.95, 20)   # Höhenlinien (Fernfeld: |psi| <= 1, im Nahfeld größer)
FARBE_PLUS, FARBE_MINUS = 'red', 'blue'   # Ladungsträger +q und -q
MS = 11                   # Größe der Ladungsträger (Punkte)
OUT = 'feldlinien_hertzscher_dipol.gif'

# ------------------------------------------------------------------ Feld auf dem Gitter
x = np.linspace(-EXTENT, EXTENT, N_GRID)
z = np.linspace(-EXTENT, EXTENT, N_GRID)
X, Z = np.meshgrid(x, z)
r = np.sqrt(X**2 + Z**2) + 1e-9                # Abstand zum Dipol (Ursprung)
sin2 = (X**2) / r**2                           # sin^2(theta) = rho^2 / r^2
P = sin2 * (1 - 1j / (k * r)) * np.exp(-1j * k * r)  # -j*rho*H_phi bis auf den Faktor k I0 d/(4 pi)

# ------------------------------------------------------------------ Bilder erzeugen
frames = []
for i in range(N_FRAMES):
    wt = 2 * np.pi * i / N_FRAMES                  # Phase omega*t
    psi = np.real(P * np.exp(1j * wt))             # Höhenlinienfunktion der momentanen E-Feldlinien

    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    ax.contour(X, Z, psi, levels=LEVELS[LEVELS > 0], colors='green', linewidths=LW, linestyles='solid')
    ax.contour(X, Z, psi, levels=LEVELS[LEVELS < 0], colors='darkblue', linewidths=LW, linestyles='solid')

    # Dipol: schwarzer Stab mit zwei gegenläufig schwingenden Ladungsträgern (Geschwindigkeit ~ I(t))
    ax.plot([0, 0], [-d_plot / 2, d_plot / 2], color='k', lw=4, solid_capstyle='butt', zorder=5)
    z_plus = d_plot / 2 * np.sin(wt)
    for zq, farbe, zeichen, zo in [(-z_plus, FARBE_MINUS, '\u2212', 6), (z_plus, FARBE_PLUS, '+', 7)]:
        ax.plot(0, zq, 'o', color=farbe, ms=MS, mec='white', mew=0.8, zorder=zo)
        ax.text(0, zq, zeichen, color='white', fontsize=MS * 0.75, fontweight='bold',
                ha='center', va='center_baseline', zorder=zo)

    ax.set_aspect('equal')
    ax.set_xlim(-EXTENT, EXTENT); ax.set_ylim(-EXTENT, EXTENT)
    ax.set_xticks(np.arange(-EXTENT, EXTENT + 0.01, 0.5)); ax.set_yticks(np.arange(-EXTENT, EXTENT + 0.01, 0.5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f')); ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel(r'x / $\lambda$'); ax.set_ylabel(r'z / $\lambda$')
    ax.set_title('Hertzscher Dipol - Feldlinien')

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
fixed = [(255, 255, 255), (0, 0, 0), (255, 0, 0), (0, 0, 255), (0, 128, 0), (0, 0, 139)]   # Farben der Grafik
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
