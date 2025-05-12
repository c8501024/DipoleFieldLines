[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.15382818)


# Dipole Field Animation Visualizations

This repository contains three Python scripts for simulating and visualizing (output as gif animations) the electromagnetic fields of a Hertzian dipole (half wave dipole currently not fully supported). It includes both 2D and 3D visualizations of electric and magnetic fields as well as the energy (Poynting vector) flow.
The scripts are inspired by the animations as done by leifiphysik.
You can easily customize the code for your own purpose. The documentation is going to be reworked in the future.

---

## 📁 Files and Outputs

| Script                        | Visualization                          | Output Format                             | Notes                                                                 |
|------------------------------|----------------------------------------|-------------------------------------------|-----------------------------------------------------------------------|
| `DIPOLANIMATION_E_ARROWS.py` | E-field lines      | **Animated GIF** E-field-2D.gif | 2D   |
| `DIPOLANIMATION_H_ARROWS.py` | H-field lines      | **Animated GIF** H-field-2D.gif | 2D                       |
| `DIPOLANIMATION_EH_3D.py`    | E- and H-field lines     | **Animated GIF** EH-Field-3D.gif          | 3D |

---

## 🧠 Physics Background

The scripts simulate electromagnetic fields of:
- A **Hertzian dipole**: ideal, infinitesimally short antenna.
- A **Finite-length dipole** (Stabdipol): two current sources.

Field computations are based on:
- Time-dependent dipole moment equations
- Electric field \( \vec{E} \), magnetic field \( \vec{H} \)
- Poynting vector \( \vec{S} = \vec{E} \times \vec{H} \)
- Energy density distribution

---

## 🚀 How to Run

### 2D Animations

```bash
python DIPOLANIMATION_E_ARROWS.py
python DIPOLANIMATION_H_ARROWS.py
```

### 3D Visualization

```bash
python DIPOLANIMATION_EH_3D.py
```

> ⚠️ This creates a directory called `DipolAnimation/` with one PNG image per time step. Currently for debugging only.

---

## 📦 Dependencies
Python version:

3.12.7

Install dependencies with:

```bash
pip install -r requirements.txt
```

### `requirements.txt` content:

```
ace_tools==0.0
ace_tools_open==0.1.0
asttokens==3.0.0
certifi==2025.1.31
charset-normalizer==3.4.1
colorama==0.4.6
contourpy==1.3.2
cycler==0.12.1
decorator==5.2.1
executing==2.2.0
fonttools==4.57.0
idna==3.10
imageio==2.37.0
ipython==9.2.0
ipython_pygments_lexers==1.1.1
itables==2.3.0
jedi==0.19.2
kiwisolver==1.4.8
matplotlib==3.10.1
matplotlib-inline==0.1.7
numpy==2.2.5
packaging==25.0
pandas==2.2.3
parso==0.8.4
pillow==11.2.1
platformdirs==4.3.7
pooch==1.8.2
prompt_toolkit==3.0.51
pure_eval==0.2.3
Pygments==2.19.1
pyparsing==3.2.3
python-dateutil==2.9.0.post0
pytz==2025.2
pyvista==0.45.0
requests==2.32.3
scipy==1.15.2
scooby==0.10.0
six==1.17.0
stack-data==0.6.3
traitlets==5.14.3
typing_extensions==4.13.2
tzdata==2025.2
urllib3==2.4.0
vtk==9.4.2
wcwidth==0.2.13
```

---

## 📝 License

This project is licensed under the **MIT License**.

---

## 🙋 Citation and Acknowledgement

If you use this software in teaching, research, or publication, please cite:

> **Dominik Mair**, *Hertzian Dipole Field Animation Tools*, University of Innsbruck, 2025.  
> DOI: [10.5281/zenodo.15382818](https://doi.org/10.5281/zenodo.15382818)

```bibtex
@misc{mair2025dipole,
  author       = {Dominik Mair},
  title        = {Hertzian Dipole Field Animation Tools},
  year         = 2025,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.15382818},
  url          = {https://doi.org/10.5281/zenodo.15382818}
}
```
---

## 🖼️ Preview

Here are some example animations of the dipole field:

**EH Field (3D)**  
<video src="EH-field-3D.mp4" autoplay loop muted playsinline width="600"></video>

**E Field (2D)**  
<video src="E-field-2D.mp4" autoplay loop muted playsinline width="600"></video>

**H Field (2D)**  
<video src="H-field-2D.mp4" autoplay loop muted playsinline width="600"></video>