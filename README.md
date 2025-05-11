# Dipole Field Animation Visualizations

This repository contains three Python scripts for simulating and visualizing the electromagnetic fields of a Hertzian or short dipole. It includes both 2D and 3D visualizations of electric and magnetic fields as well as the energy (Poynting vector) flow.

---

## 📁 Files and Outputs

| Script                        | Visualization                          | Output Format                             | Notes                                                                 |
|------------------------------|----------------------------------------|-------------------------------------------|-----------------------------------------------------------------------|
| `DIPOLANIMATION_E_ARROWS.py` | E-field lines      | **Animated GIF** (created via Matplotlib) | 2D   |
| `DIPOLANIMATION_H_ARROWS.py` | H-field lines      | **Animated GIF** (created via Matplotlib) | 2D                       |
| `DIPOLANIMATION_EH_3D.py`    | E- and H-field lines     | **Folder of rendered frames**             | Saves 3D images to `DipolAnimation/`; user must manually compile to GIF. |

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

### Run the 2D Animations:
```bash
python DIPOLANIMATION_E_ARROWS.py
python DIPOLANIMATION_H_ARROWS.py
