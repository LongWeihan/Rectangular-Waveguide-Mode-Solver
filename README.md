
# Rectangular Waveguide Mode Solver

This MATLAB-based project analyzes the guided modes in rectangular optical waveguides. It computes the effective refractive indices (Neff) for both TE and TM modes under different waveguide geometries.

## 🔍 Overview

The project contains four separate scripts, each focusing on a specific configuration:

1. **TE mode (square waveguide, a:b = 1:1)**
2. **TE mode (rectangular waveguide, a:b = 1:2)**
3. **TM mode (square waveguide, a:b = 1:1)**
4. **TM mode (rectangular waveguide, a:b = 1:2)**

Each script scans a range of waveguide dimensions and determines:
- Supported modes (based on lateral and vertical orders)
- Effective refractive indices for each mode
- Single-mode and dual-mode operation regions
- Mode dispersion curves

## 📈 Output

The scripts generate:
- Dispersion plots of effective refractive index vs waveguide size
- Highlighted regions showing single-mode and dual-mode operation
- Console output summarizing:
  - Range of single-mode operation
  - Range of dual-mode operation
  - Neff for the center of the single-mode region

## 🧠 Methods

The calculations are based on solving the characteristic equations for TE and TM modes in a rectangular waveguide using `fzero` root-finding.

- **Lateral direction (TE):** Solved using the transverse-electric wave equation
- **Vertical direction (TM):** Solved using the transverse-magnetic wave equation

## 🧰 Requirements

- MATLAB R2020+ or compatible
- No additional toolboxes required

## 📁 File Structure

```plaintext
├── te11.m        # TE mode, square waveguide (a:b = 1:1)
├── te12.m        # TE mode, rectangular waveguide (a:b = 1:2)
├── tm11.m        # TM mode, square waveguide (a:b = 1:1)
├── tm12.m        # TM mode, rectangular waveguide (a:b = 1:2)
└── README.md
```

## ✨ Example Plot
![image](https://github.com/user-attachments/assets/a3f30aa0-c615-4c94-aed5-5b9e8292bdce)


## 📜 License

This project is released under the MIT License. Feel free to use, modify, and share.

## 🤝 Contributions

Contributions and suggestions are welcome! Feel free to open an issue or pull request.

