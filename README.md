# Dipole Field Line Animations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15712770.svg)](https://doi.org/10.5281/zenodo.15712770)

Python animations of the instantaneous electric and magnetic field lines of a **Hertzian dipole** and a **half-wave dipole**. Each antenna has a 2D electric-field animation, a 2D magnetic-field animation, and a combined 3D view.

The scripts evaluate analytical field expressions and export looping GIFs. They are intended for teaching and for exploring how the near field evolves into an outward-propagating wave.

## Scripts and output files

| Script | View | Output |
|---|---|---|
| `feldlinien_hertzscher_dipol.py` | E field in the plane `y = 0` | `feldlinien_hertzscher_dipol.gif` |
| `feldlinien_halbwellendipol.py` | E field in the plane `y = 0` | `feldlinien_halbwellendipol.gif` |
| `h_feldlinien_hertzscher_dipol.py` | H field in the plane `z = 0` | `h_feldlinien_hertzscher_dipol.gif` |
| `h_feldlinien_halbwellendipol.py` | H field in the plane `z = 0` | `h_feldlinien_halbwellendipol.gif` |
| `eh_feldlinien_3d_hertzscher_dipol.py` | E and H fields in two perpendicular planes | `eh_feldlinien_3d_hertzscher_dipol.gif` |
| `eh_feldlinien_3d_halbwellendipol.py` | E and H fields in two perpendicular planes | `eh_feldlinien_3d_halbwellendipol.gif` |

`Lehrvortrag_Halbwellendipol.pptx` contains the accompanying lecture presentation. `LIESMICH.txt` provides a short German quick start.

## Installation and execution

Install the required packages in the Python environment used to run the scripts:

```powershell
python -m pip install numpy matplotlib pillow
```

The 3D scripts also import `contourpy`, which is installed as a dependency of Matplotlib. They use Matplotlib's `mplot3d` renderer. Python 3.12 is a suitable starting environment.

Open a terminal in the repository directory and run any script, for example:

```powershell
python feldlinien_halbwellendipol.py
python h_feldlinien_halbwellendipol.py
python eh_feldlinien_3d_halbwellendipol.py
```

To generate all six animations in PowerShell:

```powershell
Get-ChildItem -File -Filter '*dipol.py' | ForEach-Object {
    python $_.FullName
}
```

GIFs are written to the **current working directory**, using the filename specified by `OUT`. An existing file with that name is overwritten. The scripts use the noninteractive `Agg` backend, so no plot window opens. Rendering and palette conversion can take several minutes and use substantial memory because all frames are held in memory.

By default, each GIF contains 80 frames with 30 ms per frame: one simulated oscillation is displayed over 2.4 seconds. This is an illustrative playback speed, not the physical RF period.

## Field model and conventions

The dipole lies on the z axis in free space. Cylindrical coordinates are denoted by $(\rho,\phi,z)$, with $\rho=\sqrt{x^2+y^2}$, and $k=2\pi/\lambda$. Lengths are displayed relative to the wavelength; the scripts set `lam = 1.0`.

Complex field amplitudes use the time convention

$$
\mathbf E(\mathbf r,t)=\operatorname{Re}\{\widetilde{\mathbf E}(\mathbf r)e^{j\omega t}\},
\qquad
\mathbf H(\mathbf r,t)=\operatorname{Re}\{\widetilde{\mathbf H}(\mathbf r)e^{j\omega t}\}.
$$

Thus $e^{-jkr}$ describes outward propagation. The variable `wt` is the phase $\omega t$. The time factor is applied **once** when obtaining an instantaneous field from its complex amplitude.

The Hertzian model assumes an infinitesimally short current element with uniform current. The half-wave model assumes a thin wire of length $L=\lambda/2$ with the prescribed sinusoidal current distribution

$$
\widetilde I(z)=I_0\cos(kz),\qquad |z|\leq\lambda/4.
$$

The closed-form fields correspond to these idealized source models. Wire radius, a finite feed gap, material losses, and the actual current distribution of a particular antenna are not solved numerically.

## How the electric field lines are obtained

### 1. Use axial symmetry and Maxwell's equation

For these models, the magnetic field has only an azimuthal component $\widetilde H_\phi$, while the electric field lies in a meridional plane. Outside the source, Ampere's law is

$$
\nabla\times\widetilde{\mathbf H}=j\omega\varepsilon_0\widetilde{\mathbf E}.
$$

Writing out the curl in cylindrical coordinates gives

$$
\widetilde E_\rho=-\frac{1}{j\omega\varepsilon_0\rho}
\frac{\partial(\rho\widetilde H_\phi)}{\partial z},
\qquad
\widetilde E_z=\frac{1}{j\omega\varepsilon_0\rho}
\frac{\partial(\rho\widetilde H_\phi)}{\partial\rho}.
$$

### 2. Construct a real stream function

Define

$$
\psi(\rho,z,t)=\operatorname{Re}\{-j\rho\widetilde H_\phi(\rho,z)e^{j\omega t}\}.
$$

Then the instantaneous electric field satisfies

$$
E_\rho=-\frac{1}{\omega\varepsilon_0\rho}\frac{\partial\psi}{\partial z},
\qquad
E_z=\frac{1}{\omega\varepsilon_0\rho}\frac{\partial\psi}{\partial\rho}.
$$

Consequently,

$$
\mathbf E\cdot\nabla\psi=0.
$$

The field is tangent to curves of constant $\psi$: **contour lines of this stream function are instantaneous E-field lines**, wherever the field is regular and nonzero. This avoids integrating a separate field-line differential equation. The factor $-j$ is essential; taking contours of the instantaneous $\rho H_\phi$ alone would generally give the wrong time phase.

Real constant prefactors can be dropped when drawing the family of curves, with contour levels rescaled accordingly. A complex phase factor cannot be discarded in the same way.

### 3. Evaluate the stream function for each antenna

For the **Hertzian dipole**, with current-element length $d$, $r=\sqrt{\rho^2+z^2}$ and $\sin\theta=\rho/r$,

$$
\widetilde H_\phi=\frac{k^2 I_0d}{4\pi}
\left(\frac{j}{kr}+\frac{1}{(kr)^2}\right)e^{-jkr}\sin\theta.
$$

After dropping a real constant, its stream function becomes

$$
\psi\propto\sin^2\theta\left[
\cos(\omega t-kr)+\frac{\sin(\omega t-kr)}{kr}
\right].
$$

This is implemented by the complex array `P` in `feldlinien_hertzscher_dipol.py` and by `psi_E()` in its 3D counterpart:

```python
r = np.sqrt(X**2 + Z**2) + 1e-9
P = (X**2 / r**2) * (1 - 1j / (k * r)) * np.exp(-1j * k * r)
psi = np.real(P * np.exp(1j * wt))
```

For the **half-wave dipole**, let $a=\lambda/4$ and

$$
R_1=\sqrt{\rho^2+(z-a)^2},\qquad
R_2=\sqrt{\rho^2+(z+a)^2}.
$$

The prescribed sinusoidal current gives

$$
\widetilde H_\phi=\frac{jI_0}{4\pi\rho}
\left(e^{-jkR_1}+e^{-jkR_2}\right),
$$

so the stream function, apart from a real constant, is

$$
\psi=\cos(\omega t-kR_1)+\cos(\omega t-kR_2).
$$

In `feldlinien_halbwellendipol.py`, this is evaluated as

```python
rho = np.abs(X) + 1e-9
R1 = np.sqrt(rho**2 + (Z - h)**2)
R2 = np.sqrt(rho**2 + (Z + h)**2)
P = np.exp(-1j * k * R1) + np.exp(-1j * k * R2)
psi = np.real(P * np.exp(1j * wt))
```

Here `h = lam / 4`. Using `abs(X)` mirrors the axially symmetric solution into the left and right halves of the plane `y = 0`.

### 4. Draw contours at each time step

The 2D scripts pass `psi` to `ax.contour()` at the levels specified by `LEVELS`. Positive levels are green and negative levels are blue. These colors label the **stream-function level**, not positive and negative electric charge or the direction of the E vector.

The 3D scripts use `contourpy.contour_generator()` and `gen.lines(lev)` to obtain the same kind of contour as coordinate arrays, then place those curves in the plane `y = 0`.

For equally spaced contour levels, the local line density in the meridional drawing scales with $|\nabla\psi|=\omega\varepsilon_0\rho|\mathbf E|$. It therefore does not directly represent $|\mathbf E|$. The small `1e-9` offsets avoid division by zero; they do not model a finite wire radius or resolve source singularities.

## How the magnetic field lines are obtained

Since $\mathbf H=H_\phi\mathbf e_\phi$, its field lines at a fixed time are **circles around the z axis**. Their geometry is known directly. The scripts choose which radii to draw and determine the circulation direction from the instantaneous field.

In the equatorial plane `z = 0`, the complex amplitudes used by `H_phi()` are, up to real constants,

$$
\widetilde H_\phi(\rho)\propto
\begin{cases}
\left(\dfrac{j}{k\rho}+\dfrac{1}{(k\rho)^2}\right)e^{-jk\rho}, & \text{Hertzian dipole},\\[6pt]
\dfrac{j}{\rho}e^{-jk\sqrt{\rho^2+(\lambda/4)^2}}, & \text{half-wave dipole}.
\end{cases}
$$

For visibility over a large radial range, `h_ringe()` uses the locally normalized instantaneous field

$$
h(\rho,t)=\frac{\operatorname{Re}\{\widetilde H_\phi(\rho)e^{j\omega t}\}}{|\widetilde H_\phi(\rho)|}
=\cos\!\left(\omega t+\arg\widetilde H_\phi(\rho)\right).
$$

The implementation proceeds as follows:

1. Sample 20,000 radii and compute `h` from `H_phi()`.
2. Integrate $|h|$ radially using the trapezoidal rule to obtain a cumulative coordinate $G(\rho)$.
3. Find the zeros of `h`, dividing the radial domain into bands of constant circulation direction.
4. Place rings at equal increments in $G$, using interpolation to recover their radii. Within each band, positions are anchored to a zero of `h`; the truncated outermost band uses its inner boundary.
5. Return each radius together with the sign of `h`.

The increment is `delta = (lam / np.pi) / n_halbwelle`. In the far field, the integral of $|\cos(k\rho-\omega t)|$ over a radial half-wavelength is $\lambda/\pi$, so `n_halbwelle` controls the number of rings per half-wavelength there.

This gives a radial ring density approximately proportional to $|h|$, with the pattern evolving outward. **The radial amplitude decay of the unnormalized H field is intentionally removed.** The moving rings are a visualization of the field pattern, not material objects transported by the wave.

In 2D, triangles on the x axis indicate circulation: positive $H_\phi$ is counterclockwise when viewed from positive z toward the origin; negative $H_\phi$ is clockwise. The circles themselves are red in both cases. In the 3D view, all H circles are orange and do not have direction markers.

## Combined 3D view and charge markers

The combined view places E-field curves in `y = 0` and H-field circles in `z = 0`. It shows two representative perpendicular sections, not a volume filled with field lines. `segmente()` splits curves for manual drawing order, since Matplotlib's 3D renderer does not provide a full depth buffer. Distant H circles are faded using `H_AUSBLENDEN`.

The black rod marks the antenna. The red and blue charge markers move symbolically as

$$
z_+(t)=a_{\mathrm{plot}}\sin(\omega t),\qquad
z_-(t)=-a_{\mathrm{plot}}\sin(\omega t).
$$

Their velocities are proportional to $\cos(\omega t)$, matching the chosen current phase. Maximum separation represents maximum charge accumulation; passage through the center represents maximum current. These markers are illustrative and are not a particle simulation or the charge distribution used to compute the fields. For the Hertzian dipole, the visible rod length is also only a drawing aid.

## Customization

Edit the parameter block near the top of the relevant script:

| Parameter | Effect |
|---|---|
| `N_FRAMES` | Number of phase samples per oscillation; default 80 |
| `FRAME_MS` | GIF display time per frame; default 30 ms |
| `N_GRID` | Resolution of the E-field contour grid; 601 in 2D, 801 in 3D |
| `LEVELS`, `E_LEVELS` | E-stream-function contour levels and number of curves |
| `LINIEN_PRO_HALBWELLE` | H-ring density; default 16 in 2D, 8 in 3D |
| `EXTENT` | Display extent of the 2D views |
| `SICHTBEREICH`, `RECHENBEREICH` | 3D axis limits and larger E-field calculation domain |
| `ELEV`, `AZIM`, `ZOOM` | 3D camera settings |
| `LW`, `LW_E`, `LW_H_VORNE`, `LW_H_HINTEN` | Line widths |
| `MS`, `GROESSE_Q` | Size of the illustrative charge markers |
| `OUT` | Output GIF path |

Increasing `N_FRAMES` without changing `FRAME_MS` slows playback. To retain the same playback period, adjust both so their product stays constant, subject to GIF's 10 ms timing resolution. Higher grid resolution and more frames increase rendering time and memory use.

The export code uses a common palette across all frames, with reserved drawing colors to retain small details such as the charge markers. White margins are cropped in the 2D exports; the 3D exports retain the full canvas.

## Citation and license

Dominik Mair, *Hertzian Dipole Field Animation Tools*, 2025.

The [archived release on Zenodo](https://zenodo.org/records/15712770) has DOI [10.5281/zenodo.15712770](https://doi.org/10.5281/zenodo.15712770) and is listed under the MIT License. That archive describes an earlier release; its contents need not match the six scripts documented here. When citing these scripts, also identify the repository revision used.

```bibtex
@misc{mair2025dipole,
  author    = {Dominik Mair},
  title     = {Hertzian Dipole Field Animation Tools},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.15712770},
  url       = {https://doi.org/10.5281/zenodo.15712770}
}
```
