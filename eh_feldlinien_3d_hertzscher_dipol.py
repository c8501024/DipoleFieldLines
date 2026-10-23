# -*- coding: utf-8 -*-
"""
3D-Animation der E- und H-Feldlinien eines Hertzscher Dipols (Ansicht angelehnt an EH-field-3D aus dem
Repo DipoleFieldLines).

  * E-Feldlinien in der vertikalen Ebene y = 0 (grün/blau, wie in den 2D-Animationen):
    Höhenlinien von psi = Re{ -j * rho * H_phi * e^{j w t} } – E liegt tangential zu diesen Linien
    (der Faktor -j stammt aus rot H = j w eps0 E), mit
    -j*rho*H_phi ~ sin^2(theta) * (1 - j/(k r)) * e^{-j k r}
  * H-Feldlinien in der Äquatorebene z = 0 (orange-rote Kreise), mit
    H_phi(rho) ~ ( j/(k rho) + 1/(k rho)^2 ) * e^{-j k rho}          (theta = 90°)
    Auf die Einhüllende normiert, Dichte ~ |cos(w t + arg H_phi)|, laufen mit der Welle nach außen
    (gleiches Verfahren wie in h_feldlinien_hertzscher_dipol.py).
  * Dipol: Stab mit zwei stellvertretenden Ladungsträgern +q (rot) und -q (blau), die gegenläufig
    schwingen: z_+- = +-STAB sin(w t), Geschwindigkeit ~ I(t) = I0 cos(w t) (wie in den 2D-Animationen).

Die Felder werden über einen größeren Bereich (RECHENBEREICH) berechnet, als im Bild zu sehen ist,
damit die Feldlinien bis an den Bildrand reichen. Die H-Kreise werden nach außen hin ausgeblendet.

Umgesetzt mit matplotlib (mplot3d), damit keine weiteren Pakete nötig sind. Da matplotlib keinen
echten Tiefenpuffer hat, werden die Linien in der Reihenfolge hinten -> vorne gezeichnet.

Benötigt: numpy, matplotlib, pillow (contourpy wird mit matplotlib installiert)
"""
import io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import contourpy
from PIL import Image

# ------------------------------------------------------------------ Parameter
lam = 1.0                           # Wellenlänge (alle Längen in Vielfachen von lambda)
k = 2 * np.pi / lam                 # Wellenzahl
N_FRAMES = 80                       # Bilder pro Periode
FRAME_MS = 30                       # Anzeigedauer je Bild in ms (30 ms ~ 33 Bilder/s)
SICHTBEREICH = 2.0                  # Achsengrenzen der 3D-Ansicht (+- in lambda)
RECHENBEREICH = 4.0                 # Felder werden bis hierhin berechnet (> Sichtbereich, kein Abschneiden)
N_GRID = 801                        # Rechengitter für die E-Feldlinien
E_LEVELS = np.linspace(-0.95, 0.95, 20)   # Höhenlinien von psi
LINIEN_PRO_HALBWELLE = 8            # H-Kreise pro Halbwelle im Fernfeld
LW_E = 1.5                          # Linienstärke der E-Feldlinien
LW_H_VORNE, LW_H_HINTEN = 1.7, 1.3  # Linienstärke der H-Kreise (vordere / hintere Hälfte)
H_AUSBLENDEN = (2.4, 4.0)           # H-Kreise ab 2,4 lambda nach außen bis 4 lambda ausblenden
RHO_MIN = 0.01 * lam
STAB = 0.12 * lam                     # halbe Stablänge (nur Darstellung, gerechnet wird mit d << lambda)
GROESSE_Q = 150                     # Größe der Ladungsträger (Fläche in Punkt^2)
ELEV, AZIM, ZOOM = 20, -72, 1.9     # Blickwinkel (Grad) und Zoom
FARBE_E_POS, FARBE_E_NEG = '#008800', '#0000AA'
FARBE_H = '#FF4500'
TITEL = 'Hertzscher Dipol - E- und H-Feldlinien'
OUT = 'eh_feldlinien_3d_hertzscher_dipol.gif'


def psi_E(X, Z, wt):
    """Höhenlinienfunktion Re{-j rho H_phi e^{jwt}} in der Ebene y = 0 (Höhenlinien = momentane E-Feldlinien).

    Aus rot H = j w eps0 E folgt E ~ 1/(j w eps0 rho) * grad(rho H_phi) x e_phi; für den Momentanwert
    wandert der Faktor -j mit hinein (E läuft H im Nahfeld um eine Viertelperiode nach, im Fernfeld
    sind E und H dagegen in Phase).
    """
    r = np.sqrt(X**2 + Z**2) + 1e-9
    P = (X**2 / r**2) * (1 - 1j / (k * r)) * np.exp(-1j * k * r)
    return np.real(P * np.exp(1j * wt))


def H_phi(rho):
    """Komplexe Amplitude von H_phi in der Ebene z = 0 (ohne konstante Vorfaktoren)."""
    kr = k * rho
    return (1j / kr + 1 / kr**2) * np.exp(-1j * kr)


def h_ringe(wt, rho_max, n_halbwelle, rho_min=RHO_MIN):
    """Radien und Vorzeichen der H-Kreise zum Zeitpunkt wt.

    h = cos(w t + arg H_phi) ist das auf die Einhüllende normierte Feld. Die Kreise werden so gelegt,
    dass ihre Dichte proportional zu |h| ist; innerhalb jeder Halbwelle (zwischen zwei Nullstellen von h)
    sind sie an einer Nullstelle verankert und laufen dadurch mit der Welle nach außen.
    """
    rho = np.linspace(rho_min, rho_max, 20000)
    h = np.cos(wt + np.unwrap(np.angle(H_phi(rho))))
    G = np.concatenate([[0.0], np.cumsum(0.5 * (np.abs(h[1:]) + np.abs(h[:-1])) * np.diff(rho))])
    delta = (lam / np.pi) / n_halbwelle            # Integral von |cos| über eine Halbwelle = lambda/pi
    idx = np.where(np.sign(h[1:]) != np.sign(h[:-1]))[0]
    nodes = [rho_min] + [rho[i] - h[i] * (rho[i + 1] - rho[i]) / (h[i + 1] - h[i]) for i in idx] + [rho_max]
    ringe = []
    for a, b in zip(nodes[:-1], nodes[1:]):
        Ga, Gb = np.interp(a, rho, G), np.interp(b, rho, G)
        s = np.sign(np.interp(0.5 * (a + b), rho, h))
        if b < rho_max:      # an der äußeren Nullstelle verankern
            g = Gb - (np.arange(int((Gb - Ga) / delta + 1)) + 0.5) * delta
            g = g[g > Ga]
        else:                # äußerstes (angeschnittenes) Band an seiner inneren Nullstelle verankern
            g = Ga + (np.arange(int((Gb - Ga) / delta + 1)) + 0.5) * delta
            g = g[g < Gb]
        ringe += [(r, s) for r in np.interp(g, G, rho)]
    return ringe


def segmente(x, y, z, maske):
    """Zerlegt eine Linie in zusammenhängende Stücke, für die maske == True gilt (mit 1 Punkt Überlappung)."""
    out, start = [], None
    for i, m in enumerate(maske):
        if m and start is None:
            start = i
        if (not m or i == len(maske) - 1) and start is not None:
            end = i + 1 if m else i
            a, b = max(start - 1, 0), min(end + 1, len(maske))
            if b - a > 1:
                out.append((x[a:b], y[a:b], z[a:b]))
            start = None
    return out


def frame(wt, X, Z):
    fig = plt.figure(figsize=(10, 8), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1], projection='3d', computed_zorder=False)
    ax.view_init(elev=ELEV, azim=AZIM)
    V = SICHTBEREICH
    ax.set_xlim(-V, V); ax.set_ylim(-V, V); ax.set_zlim(-V, V)
    ax.set_box_aspect((1, 1, 1), zoom=ZOOM); ax.set_axis_off()
    # Richtung zum Betrachter: Punkte mit positivem Skalarprodukt liegen "vorne"
    v = np.array([np.cos(np.radians(ELEV)) * np.cos(np.radians(AZIM)),
                  np.cos(np.radians(ELEV)) * np.sin(np.radians(AZIM))])

    # --- E-Feldlinien (Ebene y = 0): Höhenlinien von psi, geteilt in Teile unter/über der Ebene z = 0
    gen = contourpy.contour_generator(X, Z, psi_E(X, Z, wt))
    e_unten, e_oben = [], []
    for lev in E_LEVELS:
        col = FARBE_E_POS if lev > 0 else FARBE_E_NEG
        for line in gen.lines(lev):
            x, z = line[:, 0], line[:, 1]; y = np.zeros_like(x)
            e_unten += [(s, col) for s in segmente(x, y, z, z <= 0)]
            e_oben += [(s, col) for s in segmente(x, y, z, z >= 0)]

    # --- H-Feldlinien (Ebene z = 0): Kreise, geteilt in hintere/vordere Hälfte, außen ausgeblendet
    phi = np.linspace(0, 2 * np.pi, 361)
    f0, f1 = H_AUSBLENDEN
    h_hinten, h_vorne = [], []
    for r, s in h_ringe(wt, rho_max=f1, n_halbwelle=LINIEN_PRO_HALBWELLE):
        alpha = float(np.clip((f1 - r) / (f1 - f0), 0, 1))
        if alpha <= 0:
            continue
        x, y = r * np.cos(phi), r * np.sin(phi); z = np.zeros_like(x)
        vorne = x * v[0] + y * v[1] > 0
        h_hinten += [(s_, alpha) for s_ in segmente(x, y, z, ~vorne)]
        h_vorne += [(s_, alpha) for s_ in segmente(x, y, z, vorne)]

    # --- Zeichnen von hinten nach vorne (clip_on=False: Linien reichen bis an den Bildrand)
    for (x, y, z), a in h_hinten:
        ax.plot(x, y, z, color=FARBE_H, lw=LW_H_HINTEN, alpha=0.7 * a, zorder=1, clip_on=False)
    for (x, y, z), col in e_unten:
        ax.plot(x, y, z, color=col, lw=LW_E, zorder=2, clip_on=False)
    for (x, y, z), a in h_vorne:
        ax.plot(x, y, z, color=FARBE_H, lw=LW_H_VORNE, alpha=a, zorder=3, clip_on=False)
    for (x, y, z), col in e_oben:
        ax.plot(x, y, z, color=col, lw=LW_E, zorder=4, clip_on=False)

    # --- Dipol: Stab mit zwei gegenläufig schwingenden Ladungsträgern (Geschwindigkeit ~ I(t))
    ax.plot([0, 0], [0, 0], [-STAB, STAB], color='black', lw=3.5, zorder=5, clip_on=False)
    z_plus = STAB * np.sin(wt)
    for zq, farbe, zeichen, zo in [(-z_plus, 'blue', '−', 6), (z_plus, 'red', '+', 8)]:
        ax.scatter([0], [0], [zq], color=farbe, s=GROESSE_Q, depthshade=False,
                   edgecolors='white', linewidths=0.8, zorder=zo)
        ax.text(0, 0, zq, zeichen, color='white', fontsize=9, fontweight='bold',
                ha='center', va='center', zorder=zo + 1)

    # --- Titelzeile auf weißem Balken
    fig.add_artist(Rectangle((0, 0.9), 1, 0.1, transform=fig.transFigure, color='white', zorder=10))
    fig.text(0.5, 0.955, TITEL, ha='center', va='center', fontsize=15, zorder=11)
    fig.text(0.5, 0.918, 'E: grün/blau (Ebene y = 0)    H: orange (Ebene z = 0)',
             ha='center', va='center', fontsize=10, color='#555555', zorder=11)
    buf = io.BytesIO(); fig.savefig(buf, format='png', facecolor='white'); plt.close(fig); buf.seek(0)
    return Image.open(buf).convert('RGB')


def speichere_gif(frames, out, crop=True):
    """(Optional weißen Rand abschneiden und) als GIF mit gemeinsamer Farbpalette speichern."""
    if crop:
        arrs = [np.asarray(f) for f in frames]
        mask = np.any(np.stack([(a < 245).any(axis=2) for a in arrs]), axis=0)
        ys, xs = np.where(mask)
        box = (max(xs.min() - 20, 0), max(ys.min() - 20, 0),
               min(xs.max() + 20, arrs[0].shape[1]), min(ys.max() + 20, arrs[0].shape[0]))
        frames = [f.crop(box) for f in frames]
    w, h = frames[0].size
    sheet = Image.new('RGB', (w, h * len(frames)))
    for j, f in enumerate(frames):
        sheet.paste(f, (0, j * h))
    base = sheet.resize((w // 2, h * len(frames) // 2)).quantize(colors=24, method=Image.Quantize.MEDIANCUT)
    # feste Farben der Grafik, damit seltene Farben (z. B. Ladungskugeln) nicht verloren gehen
    fixed = [(255, 255, 255), (0, 0, 0), (255, 0, 0), (0, 0, 255), (0, 136, 0), (0, 0, 170), (255, 69, 0), (153, 153, 153)]
    pal_rgb = np.array(base.getpalette()[:24 * 3] + [c for rgb in fixed for c in rgb], dtype=np.int32).reshape(-1, 3)
    pal = []
    for f in frames:   # jedes Pixel exakt der nächstgelegenen Palettenfarbe zuordnen (Weiß bleibt reinweiß)
        a = np.asarray(f).reshape(-1, 3)
        cols, inv = np.unique(a, axis=0, return_inverse=True)
        d = ((cols[:, None, :].astype(np.int32) - pal_rgb[None, :, :]) ** 2).sum(-1)
        idx = d.argmin(1).astype(np.uint8)[inv.ravel()]
        im = Image.fromarray(idx.reshape(h, w), 'P')
        im.putpalette(pal_rgb.astype(np.uint8).ravel().tolist())
        pal.append(im)
    pal[0].save(out, save_all=True, append_images=pal[1:], duration=FRAME_MS, loop=0, optimize=True)
    print('gespeichert:', out, frames[0].size)


if __name__ == '__main__':
    xs = np.linspace(-RECHENBEREICH, RECHENBEREICH, N_GRID)
    X, Z = np.meshgrid(xs, xs)
    frames = [frame(2 * np.pi * i / N_FRAMES, X, Z) for i in range(N_FRAMES)]
    speichere_gif(frames, OUT, crop=False)
