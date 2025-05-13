# Dipole Field Animation Visualizations

This repository contains three Python scripts for simulating and visualizing (output as gif animations) the electromagnetic fields of a Hertzian dipole (half wave dipole currently not fully supported). It includes both 2D and 3D visualizations of electric and magnetic fields as well as the energy (Poynting vector) flow.  For better visualization, the fields are normalized to their envelope before generating the field lines.
The scripts are inspired by the animations as done by leifiphysik.
The documentation is going to be reworked in the future. 
AI tools were used for coding.

Statement of the developer:
This code was written rather quickly with the help of AI. The main motivation behind it was the lack of easily customizable code for generating visually appealing field lines of a Hertzian dipole. Although I currently don't have the time to refactor the code for better readability, I still wanted to share it so that others can use it, create their own plots, or make use of the examples provided in the folder. Please excuse the "rough edges"—functionality came first this time!


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

Here are some example animations of the dipole field: (takes some time to load due to high quality gif)

**EH Field (3D)**  
![EH Field](Animations/EH-field-3D.gif)

**E Field (2D)**  
![E Field](Animations/E-field-2D.gif)

**H Field (2D)**  
![H Field](Animations/H-field-2D.gif)
