# VASP Band Structure Plotter

A comprehensive Python script for generating publication-quality electronic band structure plots from VASP calculations with orbital projections.
The code was applied in the peer-reviewed research article: https://doi.org/10.1103/PhysRevB.111.165137 (Figure 1.c).

## Overview

This script analyzes VASP DFT calculations and creates detailed band structure plots showing:
- Electronic band dispersions along high-symmetry k-paths
- Orbital character projections (s, p, d, f orbitals) as colored markers
- Comparison with Wannier90 interpolated bands
- Spin-resolved band structures for magnetic systems

## Features

- **Orbital Projections**: Visualize contribution of specific atomic orbitals to electronic bands
- **Multiple Plot Types**: Standard projections, VASP vs Wannier90 comparison, spin-resolved plots
- **Flexible Input**: Support for any combination of atoms and orbitals
- **Publication Ready**: High-quality PDF/PNG output with customizable formatting
- **Data Export**: Generates band_data.dat file for further analysis
- **Automatic Processing**: Handles k-path construction and Fermi level alignment

## Requirements

### Python Dependencies
```
matplotlib
numpy
termcolor
csv (built-in)
os, sys, time, math, copy (built-in)
```

### VASP Files Required
- `PROCAR` - Orbital projections
- `POSCAR` - Crystal structure
- `OUTCAR` - Calculation output
- `KPOINTS` - k-point path definition
- `INCAR` - Calculation parameters
- `DOSCAR` - Density of states
- `vasprun.xml` - Fermi level information

### Optional Files
- `wannier90_band.dat` - For Wannier90 comparison plots

## Installation

1. Ensure Python 3.x is installed with required packages:
   ```bash
   pip install matplotlib numpy termcolor
   ```

2. Place the script in your analysis directory or add to PATH

## Usage

### Basic Setup

1. **Configure the folder path** (line ~42):
   ```python
   folder = r"path/to/your/vasp/calculation"
   ```

2. **Set orbital projections** (line ~50):
   ```python
   input_print = "H.s Ni.dx2-y2"  # Hydrogen s + Nickel dx²-y² orbitals
   ```

3. **Adjust plot parameters** (lines ~20-35):
   ```python
   figure_width = 12    # cm
   figure_height = 8    # cm
   y_min = -10         # eV
   y_max = 6           # eV
   scale = 5           # marker size scaling
   ```

4. **Set Fermi level** (lines ~36-50):
   ```python
   if_manual_fermi = 1        # 1=manual, 0=read from vasprun.xml
   manual_fermi = 9.697      # eV (if manual)
   ```

### Orbital Notation

Use dot notation to specify projections:
```python
"H.s"                    # Hydrogen s orbitals
"Ni.dx2-y2"             # Nickel dx²-y² orbitals  
"O.p"                   # Oxygen p orbitals
"H.s Ni.dx2-y2"         # Multiple projections
"O.p Ni.dx2-y2 Ni.dz2"  # Up to 3 projections for RGB coloring
```

**Available orbitals:**
- **s**: s orbital
- **p**: px, py, pz (or use specific: px, py, pz)
- **d**: dxy, dyz, dz2, dxz, dx2-y2 (or use d for all)
- **f**: f orbitals (system dependent)
- **t**: total (all orbitals)

### Running the Script

```bash
python main.py
```

The script will automatically:
1. Parse all VASP files
2. Extract band structure and projections
3. Generate plots and save them as PDF/PNG
4. Export data file (if using 2 projections)

## Output Files

### Generated Plots
- **Figure_1**: Main projected band structure
- **Figure_2**: VASP vs Wannier90 comparison (if wannier90_band.dat exists)
- **Figure_3**: Spin-resolved bands (if ISPIN=2)

### File Naming Convention
```
Figure_1-[projections] h[height]w[width] [y_min]_[y_max].pdf
Example: Figure_1-H.s Ni.dx2-y2  h8w12 -10_6.pdf
```

### Data Export
- **band_data [projections].dat**: Tab-separated data file with k-points, energies, and projections

## Plot Interpretation

### Color Coding
- **2 projections**: Color scale from projection 1 (blue) → mixed (green) → projection 2 (red)
- **3 projections**: RGB coloring (Red=proj1, Green=proj2, Blue=proj3)
- **Other**: Each projection gets a different color

### Marker Size
- Larger markers = stronger orbital character
- Size controlled by `scale` parameter

### Layout
- **Black lines**: Electronic bands
- **Colored markers**: Orbital projections  
- **Vertical lines**: High-symmetry k-points
- **Horizontal line at 0**: Fermi level

## Configuration Examples

### Transition Metal Oxides
```python
input_print = "O.p Ni.dx2-y2"           # O p-bands vs Ni eg orbitals
y_min, y_max = -8, 4                     # Typical range for oxides
```

### Hydrides
```python
input_print = "H.s Ni.dx2-y2 Ni.dz2"    # H s vs Ni d orbitals
y_min, y_max = -10, 6                    # Range including H bands
```

### Heavy Elements
```python
input_print = "La.f Ni.d"               # f-d hybridization
manual_fermi = 9.67                      # Often need manual Fermi level
```

## Troubleshooting

### Common Issues

**"Error reading lmax"**
- Check OUTCAR file is complete and contains LMAX information
- Verify calculation finished successfully

**Missing k-path labels**
- Ensure KPOINTS file uses line-mode format
- Check k-path definition matches your structure

**Fermi level issues**
- Use manual Fermi level for metals: `if_manual_fermi = 1`
- Check vasprun.xml contains efermi tag

**Memory errors with large calculations**
- Reduce number of k-points or bands
- Use fewer orbital projections

**Plot formatting issues**
- Adjust figure_width, figure_height for better aspect ratio
- Modify y_min, y_max range for your system
- Change scale parameter for marker visibility

### File Requirements Check
```python
# The script checks for required files and reports missing ones
# Ensure all VASP files are in the specified folder
```

## Advanced Features

### Spin-Polarized Calculations
- Automatically handles ISPIN=2 calculations
- Creates separate spin channels if `plot_spin_channels = 1`

### Data Analysis
- Exports detailed k-point coordinates and projections
- Compatible with external analysis tools

### Custom Styling
- Modify colors array for different color schemes
- Adjust linewidth, alpha transparency
- Customize figure dimensions and DPI

## Performance Notes

- Processing time scales with number of k-points and bands
- Large projections (many atoms/orbitals) increase memory usage
- PDF generation may be slower than PNG for complex plots

## Scientific Applications

This tool is commonly used for:
- Analyzing chemical bonding in crystals
- Identifying orbital contributions to conduction/valence bands
- Comparing DFT and tight-binding models
- Studying magnetic exchange interactions
- Materials design and band gap engineering

## Citation

When using this script in publications, consider citing the relevant VASP papers and acknowledge the analysis tools used.
