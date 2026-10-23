# -*- coding: utf-8 -*-
"""
Animation der magnetischen Feldlinien (H-Feld) eines Hertzscher Dipols in der Äquatorebene z = 0 (Draufsicht).

H besitzt nur eine phi-Komponente – die Feldlinien sind Kreise um die Dipolachse.
In der Ebene z = 0 gilt (Phasor, konstante Vorfaktoren weggelassen):

    H_phi(rho) ~ ( j/(k rho) + 1/(k rho)^2 ) * e^{-j k rho}          (theta = 90°)

Darstellung (wie im Repo DipoleFieldLines): Das Feld wird auf seine Einhüllende normiert,
    h(rho, t) = Re{ H_phi e^{j w t} } / |H_phi| = cos( w t + arg H_phi(rho) ),
und die Kreise werden so gelegt, dass ihre Dichte proportional zu |h| ist. Innerhalb jeder Halbwelle sind
sie an einer Nullstelle von h verankert – so laufen sie mit der Welle nach außen, und im Nahfeld entstehen
neue Kreise am Dipol.

Farben: Kreise rot; Dreiecke auf der x-Achse zeigen die Richtung von H
(rot = H_phi > 0, gegen den Uhrzeigersinn von oben gesehen; blau = H_phi < 0, im Uhrzeigersinn).
Strom im Dipol: I(t) = I0 cos(w t).
Rechts daneben: der Dipol in Seitenansicht mit zwei stellvertretenden Ladungsträgern +q (rot) und
-q (blau), die gegenläufig schwingen: z_+- = +-STAB sin(w t), Geschwindigkeit ~ I(t) (wie in den
E-Feld-Animationen). Die gestrichelte Linie markiert die dargestellte Ebene z = 0.

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
lam = 1.0                           # Wellenlänge (alle Längen in Vielfachen von lambda)
k = 2 * np.pi / lam                 # Wellenzahl
N_FRAMES = 80                       # Bilder pro Periode
FRAME_MS = 30                       # Anzeigedauer je Bild in ms (30 ms ~ 33 Bilder/s)
EXTENT = 2.0                        # dargestellter Bereich: -2 lambda ... +2 lambda
LINIEN_PRO_HALBWELLE = 16           # Kreise pro Halbwelle im Fernfeld
LW = 1.25                           # Linienstärke der H-Feldlinien
RHO_MIN = 0.01 * lam                # innerster Radius (Dipol)
FARBE_KREIS = 'red'
FARBE_POS, FARBE_NEG = 'red', 'blue'
STAB = 0.12 * lam                     # halbe Stablänge (nur Darstellung, gerechnet wird mit d << lambda)
FARBE_PLUS, FARBE_MINUS = 'red', 'blue'   # Ladungsträger +q und -q
MS = 11                             # Größe der Ladungsträger (Punkte)
TITEL = 'Hertzscher Dipol - H-Feldlinien'
OUT = 'h_feldlinien_hertzscher_dipol.gif'


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


def dipol_seitenansicht(fig, wt):
    """Kleines Feld rechts neben dem Plot: Dipol von der Seite mit den beiden Ladungsträgern."""
    ax = fig.add_axes([0.83, 0.28, 0.14, 0.44])
    ax.set_xlim(-0.2, 0.2); ax.set_ylim(-0.42, 0.42); ax.set_aspect('equal'); ax.axis('off')
    # (gleicher Maßstab für beide Dipole: Achsen +-0.42 lambda)
    ax.plot([-0.19, 0.19], [0, 0], color='#999999', lw=0.8, ls='--', zorder=1)      # Ebene z = 0
    ax.text(0.2, 0.015, 'z = 0', color='#999999', fontsize=7, ha='right', va='bottom')
    ax.plot([0, 0], [-STAB, STAB], color='k', lw=4, solid_capstyle='butt', zorder=2)
    z_plus = STAB * np.sin(wt)
    for zq, farbe, zeichen, zo in [(-z_plus, FARBE_MINUS, '−', 3), (z_plus, FARBE_PLUS, '+', 4)]:
        ax.plot(0, zq, 'o', color=farbe, ms=MS, mec='white', mew=0.8, zorder=zo)
        ax.text(0, zq, zeichen, color='white', fontsize=MS * 0.75, fontweight='bold',
                ha='center', va='center_baseline', zorder=zo)
    ax.text(0, STAB + 0.07, 'Dipol\n(Seitenansicht)', fontsize=9, ha='center', va='bottom', clip_on=False)


if __name__ == '__main__':
    phi = np.linspace(0, 2 * np.pi, 361)
    frames = []
    for i in range(N_FRAMES):
        wt = 2 * np.pi * i / N_FRAMES                 # Phase omega*t
        fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
        ringe = h_ringe(wt, rho_max=EXTENT * np.sqrt(2), n_halbwelle=LINIEN_PRO_HALBWELLE)
        for r, s in ringe:
            ax.plot(r * np.cos(phi), r * np.sin(phi), color=FARBE_KREIS, lw=LW)
        # Richtungsmarker auf der x-Achse: rechts zeigt phi-Dach nach +y, links nach -y
        for r, s in ringe:
            if r > EXTENT:
                continue
            col = FARBE_POS if s > 0 else FARBE_NEG
            ax.plot(r, 0, '^' if s > 0 else 'v', color=col, ms=8, mec='white', mew=0.4, zorder=4)
            ax.plot(-r, 0, 'v' if s > 0 else '^', color=col, ms=8, mec='white', mew=0.4, zorder=4)
        ax.plot(0, 0, 'o', color='black', ms=10, zorder=5)            # Dipol (Draufsicht)
        ax.set_aspect('equal'); ax.set_xlim(-EXTENT, EXTENT); ax.set_ylim(-EXTENT, EXTENT)
        ax.set_xticks(np.arange(-EXTENT, EXTENT + 0.01, 0.5)); ax.set_yticks(np.arange(-EXTENT, EXTENT + 0.01, 0.5))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f')); ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.set_xlabel(r'x / $\lambda$'); ax.set_ylabel(r'y / $\lambda$'); ax.set_title(TITEL)
        dipol_seitenansicht(fig, wt)
        buf = io.BytesIO(); fig.savefig(buf, format='png', facecolor='white'); plt.close(fig); buf.seek(0)
        frames.append(Image.open(buf).convert('RGB'))
    speichere_gif(frames, OUT)
