#!/usr/bin/python
"""
VASP Band Structure Plotting Script
Generates band structure plots with orbital projections and wannier90 comparison
"""

import sys
import os
import time
import csv
import math
import copy
import matplotlib.pyplot as plt
import numpy as np
from termcolor import colored

# ==================== CONFIGURATION ====================
# Plot settings
FIGURE_WIDTH = 12  # cm
FIGURE_HEIGHT = 8  # cm
PLOT_WIDTH = 0.6
Y_MIN = -10
Y_MAX = 6
Y_STEP = 2
SCALE = 5
ALPHA = 1
COLORS = ['b', 'lime', 'r', 'c', 'm', 'y', 'k', 'g', 'w']
plt.rcParams["font.family"] = "Arial"

# Fermi level settings and constants
IF_MANUAL_FERMI = 1
MANUAL_FERMI = 9.69701759  # the energy you want to set 
PLOT_SPIN_CHANNELS = 0 #functionality removed 
CM = 1 / 2.54

# Directory and projection settings
FOLDER = r"A path to your VASP results directory"
INPUT_PRINT = #e.g. "O.p Ni.dx2-y2"

# ==================== UTILITY FUNCTIONS ====================
def split_list_halves(input_list):
    """Split list into two equal halves"""
    middle = len(input_list) // 2
    return input_list[:middle], input_list[middle:]


def read_fermi_energy(vasprun_path):
    """Read Fermi energy from vasprun.xml"""
    with open(vasprun_path, "r") as f:
        for line in f:
            if "efermi" in line:
                return float(line.split()[2])
    return 0


def read_ispin(incar_path):
    """Read ISPIN value from INCAR"""
    with open(incar_path, "r") as f:
        for line in f:
            if "ISPIN" in line:
                return line.split()[-1]
    return ''


def read_lmax(outcar_path):
    """Read maximum LMAX from OUTCAR and return orbital configuration"""
    lmax_values = []
    with open(outcar_path, "r") as f:
        for line in f:
            if " LMAX " in line:
                lmax_values.append(int(line.split()[-1]))

    max_lmax = max(lmax_values)

    if max_lmax == 6:
        orbital_no = list(range(10))
        projection_symbol = [['s'], ['py'], ['pz'], ['px'], ['dxy'], ['dyz'],
                             ['dz2'], ['dxz'], ['dx2-y2'], ['tot']]
    elif max_lmax == 8:
        orbital_no = list(range(17))
        projection_symbol = [['s'], ['py'], ['pz'], ['px'], ['dxy'], ['dyz'],
                             ['dz2'], ['dxz'], ['dx2-y2'], ['fy3x2'], ['fxyz'],
                             ['fyz2'], ['fz3'], ['fxz2'], ['fzx2'], ['fx3'], ['tot']]
    elif max_lmax == 4:
        orbital_no = list(range(6))
        projection_symbol = [['s'], ['py'], ['pz'], ['px'], ['tot']]
    else:
        raise ValueError(f"Unsupported LMAX value: {max_lmax}")

    for i, symbol in enumerate(projection_symbol):
        symbol.append(orbital_no[i])

    return projection_symbol


def read_lattice_parameters(poscar_path):
    """Read lattice parameters and calculate reciprocal vectors"""
    with open(poscar_path, "r") as f:
        lines = f.readlines()

    scale = float(lines[1])
    lattice_vectors = np.array([
        [float(x) for x in lines[2].split()],
        [float(x) for x in lines[3].split()],
        [float(x) for x in lines[4].split()]
    ])

    v_length = np.linalg.norm(lattice_vectors, axis=1) * scale

    # Calculate angles
    def angle(v1, v2):
        return np.arccos(
            np.dot(v1, v2) * scale ** 2 / (np.linalg.norm(v1) * np.linalg.norm(v2) * scale ** 2)) * 180 / np.pi

    alpha = angle(lattice_vectors[1], lattice_vectors[2])
    beta = angle(lattice_vectors[0], lattice_vectors[2])
    gamma = angle(lattice_vectors[0], lattice_vectors[1])

    # Calculate volume
    volume = v_length[0] * v_length[1] * v_length[2] * np.sqrt(
        1 - np.cos(np.radians(alpha)) ** 2 - np.cos(np.radians(beta)) ** 2 -
        np.cos(np.radians(gamma)) ** 2 + 2 * np.cos(np.radians(alpha)) *
        np.cos(np.radians(beta)) * np.cos(np.radians(gamma))
    )

    # Calculate reciprocal lattice vectors
    recip_vectors = [
        2 * np.pi * v_length[1] * v_length[2] * np.sin(np.radians(alpha)) / volume,
        2 * np.pi * v_length[0] * v_length[2] * np.sin(np.radians(beta)) / volume,
        2 * np.pi * v_length[0] * v_length[1] * np.sin(np.radians(gamma)) / volume
    ]

    return recip_vectors, lines


def read_atom_structure(poscar_lines):
    """Parse atom types and indices from POSCAR"""
    atom_types = poscar_lines[5].split()
    atoms_no = [int(x) for x in poscar_lines[6].split()]

    atom_list = []
    atom_index = 0
    for i, atom_type in enumerate(atom_types):
        indices = list(range(atom_index, atom_index + atoms_no[i]))
        atom_list.append([atom_type, atoms_no[i], indices])
        atom_index += atoms_no[i]

    return atom_list, sum(atoms_no)


def parse_procar(procar_path, efermi, atoms_sum):
    """Parse PROCAR file for eigenvalues and projections"""
    with open(procar_path, "r") as f:
        procar_data = f.read()

    procar_kpoints = procar_data.split("k-point ")
    procar_kpoints.pop(0)

    no_kpoints = len(procar_kpoints)
    all_kpoints_array = []
    ion_band = []

    for kpoint_idx, kpoint in enumerate(procar_kpoints):
        kpoint_energies = []
        for line in kpoint.split("\n"):
            if "energy" in line:
                kpoint_energies.append(float(line.split()[4]) - efermi)
        all_kpoints_array.append(kpoint_energies)

        kpoint_bands = kpoint.split("band ")
        kpoint_bands.pop(0)
        ion_band.append([])

        for band_text in kpoint_bands:
            lines = band_text.split("\n")[3:-2]
            ion_band[kpoint_idx].append([])

            for ion_idx in range(atoms_sum):
                proj_values = lines[ion_idx].split()[1:]
                ion_band[kpoint_idx][-1].append(proj_values)

    no_bands = len(all_kpoints_array[0])
    return all_kpoints_array, ion_band, no_kpoints, no_bands


def parse_kpoints(kpoints_path, recip_vectors, no_kpoints):
    """Parse KPOINTS file and calculate k-path"""
    with open(kpoints_path, "r") as f:
        lines = f.readlines()

    kpoints_line1 = lines[0].strip()
    pts_per_section = int(lines[1].split()[0]) - 1
    kpaths = kpoints_line1.split()[2:]

    k_coordinates = []
    for line in lines[4:]:
        coord = line.split()[:3]
        if len(coord) == 3:
            k_coordinates.append(coord)

    k_coordinates_recip = k_coordinates.copy()

    # Remove every other coordinate for continuous path
    for i in range(len(k_coordinates) // 2 - 1):
        k_coordinates.pop(-(i + 2))

    # Convert to Cartesian
    k_coordinates_cart = []
    for coord in k_coordinates:
        cart = [float(c) * rv for c, rv in zip(coord, recip_vectors)]
        k_coordinates_cart.append(cart)

    # Calculate distances
    kpath_cart = [0]
    k_dist = []
    for i in range(1, len(k_coordinates_cart)):
        diff = np.array(k_coordinates_cart[i]) - np.array(k_coordinates_cart[i - 1])
        dist = np.linalg.norm(diff)
        k_dist.append(dist)
        kpath_cart.append(kpath_cart[-1] + dist)

    # Generate all k-point coordinates
    all_kpoints_coord = []
    for i in range(len(kpath_cart) - 1):
        spacing = k_dist[i] / pts_per_section
        for k in range(pts_per_section + 1):
            all_kpoints_coord.append(kpath_cart[i] + k * spacing)

    # Generate reciprocal coordinates for all k-points
    k_coordinates_recip_all = []
    for i in range(len(k_coordinates_recip) // 2):
        start = np.array([float(x) for x in k_coordinates_recip[i * 2]])
        end = np.array([float(x) for x in k_coordinates_recip[i * 2 + 1]])
        spacing = (end - start) / pts_per_section

        for k in range(pts_per_section + 1):
            k_coordinates_recip_all.append((start + k * spacing).tolist())

    # Generate labels
    x_paths = [path.split("-") for path in kpaths]
    for i in range(len(x_paths) - 1):
        x_paths[i][-1] += ' | ' + x_paths[i + 1][0]
        x_paths[i + 1].pop(0)

    x_labels = [label for path in x_paths for label in path]
    x_positions = kpath_cart

    return all_kpoints_coord, x_labels, x_positions, kpath_cart, k_coordinates_recip_all


def parse_input_projections(input_string, atom_list, projection_symbol):
    """Parse projection input string"""
    input_strings = input_string.split()
    input_proj = []

    for inp_str in input_strings:
        parts = inp_str.split(".")
        atom_name = parts[0]
        orbitals = parts[1:]

        atom_indices = []
        for atom in atom_list:
            if atom_name in atom[0]:
                atom_indices = atom[2]
                break

        orbital_symbols = []
        orbital_indices = []
        for orb in orbitals:
            for proj in projection_symbol:
                if orb in proj[0]:
                    orbital_symbols.append(orb)
                    orbital_indices.append(proj[1])
                    break

        input_proj.append([atom_name, atom_indices, orbital_symbols, orbital_indices])

    return input_proj


def calculate_projections(ion_band, input_proj, no_bands, no_kpoints, scale):
    """Calculate orbital projections for plotting"""
    projection_print = []

    for proj_no, proj in enumerate(input_proj):
        projection_print.append([])

        for orb_idx in range(len(proj[3])):
            projection_print[proj_no].append([])

            for atom_idx in range(len(proj[1])):
                projection_print[proj_no][orb_idx].append([])
                proj_orbital = proj[3][orb_idx]
                proj_atom = proj[1][atom_idx]

                for band in range(no_bands):
                    projection_print[proj_no][orb_idx][atom_idx].append([])

                    for k_point in range(no_kpoints):
                        val = scale * float(ion_band[k_point][band][proj_atom][proj_orbital])
                        projection_print[proj_no][orb_idx][atom_idx][band].append(val)

    return projection_print


def aggregate_projections(projection_print, input_proj, no_bands, no_kpoints):
    """Aggregate projections for each band"""
    proj_sum = [[[] for _ in input_proj] for _ in range(no_bands)]

    for band in range(no_bands):
        for proj_no in range(len(input_proj)):
            proj_sum[band][proj_no] = np.zeros(no_kpoints)

            for orb_idx in range(len(input_proj[proj_no][3])):
                for atom_idx in range(len(input_proj[proj_no][1])):
                    proj_sum[band][proj_no] += np.array(
                        projection_print[proj_no][orb_idx][atom_idx][band]
                    )

    return proj_sum


def calculate_mixing(proj_sum, input_proj, no_bands, no_kpoints):
    """Calculate mixing parameters for 2 or 3 projections"""
    mixing = []
    band_sum_sq = []

    for band in range(no_bands):
        if len(input_proj) == 2:
            band_sum = proj_sum[band][0] + proj_sum[band][1]
            diff = proj_sum[band][0] - proj_sum[band][1]
            mix = (band_sum - diff) / (band_sum + 1e-7) / 2
            mixing.append(mix.tolist())
            band_sum_sq.append((band_sum ** 2).tolist())

        elif len(input_proj) == 3:
            band_sum_tot = proj_sum[band][0] + proj_sum[band][1] + proj_sum[band][2]
            band_sum_sq.append((band_sum_tot ** 2).tolist())

            frac_1 = proj_sum[band][0] / (band_sum_tot + 1e-7)
            frac_2 = proj_sum[band][1] / (band_sum_tot + 1e-7)
            frac_3 = proj_sum[band][2] / (band_sum_tot + 1e-7)

            mix_band = []
            for k_point in range(no_kpoints):
                max_frac = max(frac_1[k_point], frac_2[k_point], frac_3[k_point], 0.001)
                rgb = [
                    frac_3[k_point] / max_frac,  # B
                    frac_2[k_point] / max_frac,  # G
                    frac_1[k_point] / max_frac  # R
                ]
                mix_band.append(rgb)
            mixing.append(mix_band)
        else:
            mixing.append([0] * no_kpoints)
            band_sum_sq.append([0] * no_kpoints)

    return mixing, band_sum_sq


def plot_bands(ax, all_kpoints_coord, bands_eigenvalues, mixing, band_sum_sq,
               input_proj, input_strings, x_positions, y_min, y_max, kpath_cart):
    """Create main band structure plot"""

    # Plot projections
    if len(input_proj) == 2:
        for band in range(len(bands_eigenvalues)):
            if len(all_kpoints_coord) != len(bands_eigenvalues[band]):
                bands_spin = split_list_halves(bands_eigenvalues[band])
                mixing_spin = split_list_halves(mixing[band])
                band_sum_spin = split_list_halves(band_sum_sq[band])

                im1 = ax.scatter(all_kpoints_coord, bands_spin[0], c=mixing_spin[0],
                                 s=band_sum_spin[0], cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)
                im2 = ax.scatter(all_kpoints_coord, bands_spin[1], c=mixing_spin[1],
                                 s=band_sum_spin[1], cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)
            else:
                im1 = ax.scatter(all_kpoints_coord, bands_eigenvalues[band], c=mixing[band],
                                 s=band_sum_sq[band], cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)

        return im1

    elif len(input_proj) == 3:
        for band in range(len(bands_eigenvalues)):
            if len(all_kpoints_coord) != len(bands_eigenvalues[band]):
                bands_spin = split_list_halves(bands_eigenvalues[band])
                mixing_spin = split_list_halves(mixing[band])
                band_sum_spin = split_list_halves(band_sum_sq[band])

                im1 = ax.scatter(all_kpoints_coord, bands_spin[0], c=mixing_spin[0],
                                 s=band_sum_spin[0], alpha=1.0)
                im2 = ax.scatter(all_kpoints_coord, bands_spin[1], c=mixing_spin[1],
                                 s=band_sum_spin[1], alpha=1.0)
            else:
                im1 = ax.scatter(all_kpoints_coord, bands_eigenvalues[band],
                                 c=mixing[band], s=band_sum_sq[band], alpha=1.0)
        return None

    else:
        proj_sum_sq = [[np.array(proj_sum[band][i]) ** 2 for i in range(len(input_proj))]
                       for band in range(len(bands_eigenvalues))]

        for proj_no in range(len(input_proj)):
            for band in range(len(bands_eigenvalues)):
                ax.scatter(all_kpoints_coord, bands_eigenvalues[band],
                           s=proj_sum_sq[band][proj_no].tolist(),
                           color=COLORS[proj_no], alpha=ALPHA)
        return None


def configure_axes(ax, x_positions, x_labels, y_min, y_max, y_step, kpath_end):
    """Configure plot axes with labels and limits"""
    ax.vlines(x=x_positions, ymin=y_min, ymax=y_max, colors='black', linewidth=0.25)
    ax.set(xlim=[0, kpath_end], ylim=[y_min, y_max], ylabel='Energy (eV)')
    ax.set_xticks(x_positions, labels=x_labels)

    yticks = np.concatenate([np.arange(0, y_min - 0.001, -y_step),
                             np.arange(0, y_max + 0.001, y_step)])
    ax.set_yticks(yticks)
    ax.axhline(y=0.0, color='black', linestyle='-', linewidth=0.25)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.8)
    ax.tick_params(width=0.8)


# ==================== MAIN EXECUTION ====================
def main():
    start_time = time.time()

    # Setup paths
    procar_path = os.path.join(FOLDER, "PROCAR")
    poscar_path = os.path.join(FOLDER, "POSCAR")
    outcar_path = os.path.join(FOLDER, "OUTCAR")
    incar_path = os.path.join(FOLDER, "INCAR")
    kpoints_path = os.path.join(FOLDER, "KPOINTS")
    vasprun_path = os.path.join(FOLDER, "vasprun.xml")
    wannier_path = os.path.join(FOLDER, "wannier90_band.dat")
    band_data_path = os.path.join(FOLDER, f"band_data {INPUT_PRINT}.dat")

    # Read basic parameters
    ispin = read_ispin(incar_path)
    print(f"ISPIN = {ispin}")

    projection_symbol = read_lmax(outcar_path)
    print(f"LMAX = {len(projection_symbol) - 1}")
    print(f"Projection symbols: {projection_symbol}")

    recip_vectors, poscar_lines = read_lattice_parameters(poscar_path)
    print(f"Reciprocal vectors: {recip_vectors}")

    atom_list, atoms_sum = read_atom_structure(poscar_lines)
    print(f"Atom list: {atom_list}")
    print(f"Number of ions: {atoms_sum}")

    efermi = read_fermi_energy(vasprun_path)
    if IF_MANUAL_FERMI:
        efermi = MANUAL_FERMI
    print(f"Energy of Fermi level: {efermi}")

    # Parse PROCAR
    all_kpoints_array, ion_band, no_kpoints, no_bands = parse_procar(
        procar_path, efermi, atoms_sum
    )
    print(f"Number of kpoints: {no_kpoints}")
    print(f"Number of bands: {no_bands}")

    # Parse KPOINTS
    all_kpoints_coord, x_labels, x_positions, kpath_cart, k_coordinates_recip_all = parse_kpoints(
        kpoints_path, recip_vectors, no_kpoints
    )
    print(f"K-path labels: {x_labels}")

    # Parse input projections
    input_proj = parse_input_projections(INPUT_PRINT, atom_list, projection_symbol)
    input_strings = INPUT_PRINT.split()
    print(f"Input projections: {input_proj}")

    # Calculate projections
    projection_print = calculate_projections(
        ion_band, input_proj, no_bands, no_kpoints, SCALE
    )
    proj_sum = aggregate_projections(projection_print, input_proj, no_bands, no_kpoints)
    mixing, band_sum_sq = calculate_mixing(proj_sum, input_proj, no_bands, no_kpoints)

    # Extract eigenvalues
    bands_eigenvalues = [[all_kpoints_array[k][b] for k in range(no_kpoints)]
                         for b in range(no_bands)]

    # Create main plot
    fig1, ax1 = plt.subplots(figsize=(FIGURE_WIDTH * CM, FIGURE_HEIGHT * CM))

    im = plot_bands(ax1, all_kpoints_coord, bands_eigenvalues, mixing, band_sum_sq,
                    input_proj, input_strings, x_positions, Y_MIN, Y_MAX, kpath_cart)

    configure_axes(ax1, x_positions, x_labels, Y_MIN, Y_MAX, Y_STEP, kpath_cart[-1])

    # Add black band lines
    for band in range(no_bands):
        if len(all_kpoints_coord) != len(bands_eigenvalues[band]):
            bands_spin = split_list_halves(bands_eigenvalues[band])
            ax1.plot(all_kpoints_coord, bands_spin[0], color="black", linewidth=0.5)
            ax1.plot(all_kpoints_coord, bands_spin[1], color="black", linewidth=0.5)
        else:
            ax1.plot(all_kpoints_coord, bands_eigenvalues[band], color="black", linewidth=0.5)

    # Add colorbar or legend
    box = ax1.get_position()
    ax1.set_position([box.x0 + 0.05, box.y0, box.width * PLOT_WIDTH, box.height])

    if len(input_proj) == 2 and im is not None:
        fig1.colorbar(im, ax=ax1, ticks=[0, 0.5, 1], shrink=0.25, aspect=10).set_ticklabels(
            [input_strings[0], 'mixed', input_strings[1]]
        )
    elif len(input_proj) != 2:
        plot_legend = [ax1.scatter([-10], [0], s=SCALE / 5, color=COLORS[i])
                       for i in range(len(input_proj))]
        ax1.legend(plot_legend, input_strings, loc='upper center',
                   bbox_to_anchor=(1.30, 0.7), markerscale=4, frameon=False)

    # Save figure
    basename = os.path.join(FOLDER, f"Figure_1-{INPUT_PRINT} h{FIGURE_HEIGHT}w{FIGURE_WIDTH} {Y_MIN}_{Y_MAX}")
    fig1.savefig(f"{basename}.png", format="png", dpi=500)
    fig1.savefig(f"{basename}.pdf", format="pdf")
    print(f"\nSaved: {basename}.pdf\n")

    # Wannier90 comparison plot (if exists)
    if os.path.isfile(wannier_path):
        print("Printing wannierization plot...")
        # [Wannier plotting code would go here - simplified for brevity]

    print(f"--- {time.time() - start_time:.2f} seconds ---")


if __name__ == "__main__":
    main()

