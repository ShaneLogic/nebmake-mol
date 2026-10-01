# nebmake-mol

**Generate NEB starting structures for molecular reorientation in hybrid perovskites.**

`nebmake-mol` is a Python tool that builds a sequence of VASP structures between two endpoint configurations. It separates molecular translation and rotation from the motion of the inorganic framework, retaining the initial molecule's internal geometry during interpolation in a fixed cell.

The current implementation is configured for methylammonium lead iodide, CH3NH3PbI3 (MAPbI3), with a Pb-I framework and C-N-H molecular units.

> This tool generates an initial path. It does not run DFT, optimize a nudged elastic band (NEB), locate a transition state, or calculate an activation barrier. Inspect the generated structures before using them in a separate NEB calculation.

- [Method](#method)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input requirements](#input-requirements)
- [Output files](#output-files)
- [Verification](#verification)
- [Limitations and troubleshooting](#limitations-and-troubleshooting)
- [Code map](#code-map)

## Method

### Why separate molecular rotation from atomic displacement?

A conventional linear interpolation moves each atom independently between its endpoint positions:

$$
\mathbf{r}_i(s)=(1-s)\mathbf{r}_i^{(0)}+s\mathbf{r}_i^{(1)},\qquad 0\leq s\leq 1.
$$

For a rotating molecule, these straight trajectories generally do not preserve bond lengths or angles. For example, if a bond vector reverses direction between the endpoints, its linearly interpolated midpoint is zero. This is an unsuitable starting geometry for an intact molecular ion.

`nebmake-mol` instead describes each molecule using a reference center, an orientation, and coordinates in a frame attached to the molecule. The center and orientation change along the path; the internal coordinates are taken from the initial structure.

### 1. Identify molecules across periodic boundaries

The code selects reference atoms from the least abundant organic element in the input and collects nearby C, N, and H atoms within a user-specified cutoff. Neighbor coordinates are shifted by periodic lattice translations so that each molecular group can be represented together across a cell boundary.

For the expected `C N H` ordering and MA stoichiometry, C supplies the reference atoms. The cutoff is a **reference-atom neighbor-search radius in angstroms**, not an interpolation step or a target bond length. It must include every atom in one molecule without collecting atoms from another.

The extracted coordinates and a mass-weighted reference center are saved for each molecule. The legacy center calculation counts the reference atom twice; the resulting point should be interpreted as the implementation's reference center, rather than the conventional center of mass. The coordinate decomposition below remains valid about that reference point.

### 2. Encode translation, orientation, and internal geometry

Let the molecular reference center be $\mathbf{c}$ and the orthonormal molecular frame be

$$
R=[\mathbf{g}_1\ \mathbf{g}_2\ \mathbf{g}_3].
$$

The primary axis $\mathbf{g}_1$ points from the center toward one selected atom. A second reference atom defines a perpendicular direction $\mathbf{g}_2$, and $\mathbf{g}_3=\mathbf{g}_1\times\mathbf{g}_2$. With the expected MA ordering, these references are the N atom and a selected H atom.

Each atom is encoded by its coordinates in this molecular frame:

$$
\mathbf{q}_i=R^\mathsf{T}(\mathbf{r}_i-\mathbf{c}),
\qquad
\mathbf{r}_i=\mathbf{c}+R\mathbf{q}_i.
$$

The implementation stores the components of $\mathbf{q}_i$ as `l1`, `l2`, and `l3`. Three angles describe the orientation:

| Angle | Meaning |
| --- | --- |
| `theta` | Azimuth of the primary axis in the laboratory xy plane. |
| `phi` | Polar angle of the primary axis measured from the laboratory z axis. |
| `gamma` | Rotation of the secondary direction about the primary axis, relative to a constructed perpendicular basis. |

These are the angle conventions in `src/coordinates_convert.py` and `src/molecule_code.py`; they should not be assumed to match another package's Euler-angle convention.

### 3. Interpolate the center and orientation

For `N` output structures, the interpolation parameter is

$$
s_k=\frac{k}{N-1},\qquad k=0,\ldots,N-1.
$$

The molecular center and the selected angular coordinates are interpolated linearly. Each molecular image is then reconstructed using the **initial** internal coordinates:

$$
\mathbf{c}(s)=(1-s)\mathbf{c}_0+s\mathbf{c}_1,
\qquad
\mathbf{r}_i(s)=\mathbf{c}(s)+R\bigl(\theta(s),\phi(s),\gamma(s)\bigr)\mathbf{q}_i^{(0)}.
$$

An ideal orthonormal rotation preserves all intramolecular pair distances when the cell is fixed. In practice, the code rounds several intermediate quantities and writes coordinates with finite precision, so small numerical deviations remain.

The current implementation applies branch adjustments to `theta` and `phi`, then interpolates `theta`, `phi`, and `gamma` independently. It does not use quaternion SLERP or guarantee the shortest rotation in three dimensions. In particular, `gamma` has no analogous periodic branch adjustment. Large rotations and angular branch crossings need visual inspection.

**Internal deformation is not interpolated.** If the relaxed endpoint molecules have different bond lengths or bond angles, the reconstructed final molecular geometry can differ from the supplied final structure.

### 4. Interpolate the framework and assemble the images

The Pb and I atoms are interpolated in fractional coordinates using corresponding site indices and component-wise periodic displacement wrapping:

$$
\Delta\mathbf{f}=\mathbf{f}_1-\mathbf{f}_0-\mathrm{round}(\mathbf{f}_1-\mathbf{f}_0).
$$

$$
\mathbf{f}(s)=\mathbf{f}_0+s\Delta\mathbf{f}.
$$

Here, $\mathbf{f}_0$ and $\mathbf{f}_1$ are the endpoint fractional coordinates. The function `round` selects the nearest integer separately for each component, following NumPy's rounding convention; subtracting it wraps each fractional displacement into the interval from -0.5 to 0.5.

The lattice is interpolated using the stretch part of a polar decomposition. If $A_0$ and $A_1$ contain the lattice vectors as rows, the code forms

$$
F=A_1^\mathsf{T}(A_0^\mathsf{T})^{-1}=UP.
$$

$$
A(s)^\mathsf{T}=[I+s(P-I)]A_0^\mathsf{T}.
$$

Only the stretch $P$ is used; the rotational factor $U$ is discarded. Endpoint cells therefore need a consistent orientation and lattice-vector correspondence.

Finally, the interpolated Pb-I coordinates are combined with the reconstructed molecular coordinates. For changing cells, the molecular branch uses the initial lattice for its coordinate conversions, while the output uses the interpolated lattice. This adds an affine deformation: **strict Cartesian rigid-molecule preservation applies to the fixed-cell case**, not generally to a variable-cell path.

## Installation

The program runs directly from the source tree. There is no package installation or installed `nebmake-mol` command.

The runtime dependencies are NumPy, pandas, SciPy, and pymatgen. The following environment was used to smoke-test the documented command with a synthetic fixed-cell MA rotation:

```bash
conda create -n nebmake-mol python=3.13 pip
conda activate nebmake-mol
python -m pip install numpy==2.3.4 pandas==2.2.3 scipy==1.16.3 pymatgen==2025.10.7
```

That earlier smoke test used Python 3.13.5 and produced six readable 12-atom structures. It checks the command and file workflow, not the physical quality of an NEB path or compatibility with every input system. A fresh source-review check is described under [Verification](#verification).

The original README records pymatgen `2023.8.10`. The repository also includes `requirements-pip.txt` and `requirements-conda.txt` as historical environment records; review their version and platform constraints before using them to reconstruct that environment.

## Quick start

### 1. Create a dedicated run directory

The current code defines `base_dir` as the **parent directory of the source checkout**, regardless of the shell's working directory. Keep one checkout and one pair of endpoint copies in each run directory:

```bash
mkdir nebmake-run
cd nebmake-run
git clone https://github.com/ShaneLogic/nebmake-mol.git
```

Prepare the endpoint files according to [Input requirements](#input-requirements), then copy them into this directory. Replace the example source paths below with the paths to your prepared structures:

```bash
cp /absolute/path/to/initial.vasp initial.vasp
cp /absolute/path/to/final.vasp final.vasp
```

The layout before interpolation is:

```text
nebmake-run/
  initial.vasp
  final.vasp
  nebmake-mol/
    config.py
    interp_main.py
    src/
    utils/
```

> Run on copies: the first processing step rewrites `initial.vasp` in place to adjust periodic displacements. Use a fresh run directory for a new endpoint pair or cutoff, because existing extracted-molecule files are reused.

### 2. Generate the path

From `nebmake-run/`, with the environment activated:

```bash
python nebmake-mol/interp_main.py initial.vasp final.vasp 6 2.5
```

| Positional argument | Example | Meaning |
| --- | --- | --- |
| Initial structure | `initial.vasp` | VASP-format initial endpoint. Relative paths are resolved against `base_dir`. This file is modified in place. |
| Final structure | `final.vasp` | VASP-format final endpoint, with matching composition and atom correspondence. |
| Number of structures | `6` | Total number of generated structures, **including both endpoints**. Use an integer of at least 2. |
| Molecular cutoff | `2.5` | Radius in angstroms used to collect organic atoms around each reference atom. Check it for your geometry. |

This command produces six numbered directories, `00` through `05`. There are **four intermediate images**, so a conventional fixed-cell VASP NEB calculation using this set would use `IMAGES = 4`, after validating the endpoints and path.

Absolute input paths are also accepted, but output still goes into `base_dir`. Running the script from a different working directory does not redirect its output.

### 3. Inspect and prepare the NEB calculation

1. Visualize the complete path and confirm the intended molecular rotation and translation.
2. Check molecular membership, bond geometry, and short contacts with the framework and other molecules.
3. Compare `00/POSCAR` and the last numbered `POSCAR` with the intended endpoints, allowing for periodic translations and the generated atom order. They are reconstructed structures, not untouched copies of the input files. If exact relaxed endpoints are required, restore them only after matching every atomic index to the path.
4. Prepare the electronic-structure input files separately. Match `POTCAR` and all species-indexed settings to the generated `Pb I C N H` order.
5. Use the appropriate NEB implementation for your problem. A variable-cell path needs a solver that supports cell degrees of freedom; generating different cells does not make an ordinary fixed-cell NEB calculation a solid-state NEB calculation.

## Input requirements

The parsers and assembly routines assume a specific POSCAR layout. Use:

- VASP 5-style files with an explicit element-symbol line and element-count line.
- `Direct` fractional coordinates with exactly three numeric values per atom.
- No `Selective dynamics` line, selective-dynamics flags, or trailing atom labels in the coordinate block.
- The same atom counts, element order, and consistent atom-to-atom and molecule-to-molecule correspondence at both endpoints. The command does not solve an atom-assignment problem.
- The element order **`Pb I C N H`** at both endpoints for the default configuration. In particular, keep the organic order `C N H` consistent with the hard-coded output header.
- Consistent cell orientation, lattice-vector correspondence, and periodic images. A fixed cell is the simplest starting point.

Each extracted MA group should contain one C, one N, and six H atoms. Check the `center_neighbor/*.dat` files: the first row is the reference center, followed by the molecular atomic coordinates, all in fractional coordinates. An MA group therefore has nine rows in total.

The organic masses and species lists are defined in `config.py`. Supporting a different composition requires reviewing molecule identification, reference-atom selection, element grouping, and POSCAR assembly in addition to changing those lists. The present command is not a general molecular-connectivity detector.

## Output files

For `N = 6`, output is written next to the source checkout:

```text
nebmake-run/
  initial.vasp                 # Periodic-image-adjusted working copy
  final.vasp
  nebmake-mol/
  center_neighbor/
    initial_0.dat             # Center and coordinates of molecule 0
    final_0.dat
    ...                       # More files for additional molecules
  interp_result/
    interp_1.vasp
    ...
    interp_6.vasp
  00/POSCAR
  01/POSCAR
  02/POSCAR
  03/POSCAR
  04/POSCAR
  05/POSCAR
```

| Output | Purpose and current behavior |
| --- | --- |
| `center_neighbor/*.dat` | Extracted molecular groups. Existing files are not overwritten; stale files can affect a later run. |
| Numbered `*/POSCAR` files | Assembled structures copied into NEB-style directories before the final alignment step. Use these as the starting path for the inspection workflow above. Existing files with the same names are overwritten. |
| `interp_result/interp_*.vasp` | Assembled structures subsequently shifted to align their arithmetic mean fractional coordinates with the first structure. The shift is applied to every atom in each image. |

**The two structure sets are not necessarily identical.** `utils/reset_lattice.py` runs after the numbered directories have been populated and only changes `interp_result/`. Despite its name, this step translates coordinates; it does not reset the lattice vectors. Do not mix the two sets without checking the applied shifts.

Directory names are constructed by prefixing the integer index with `0`: index 9 is `09`, but index 10 is `010`. For more than ten total structures, rename the directories as needed for the target solver's naming convention and verify their numerical order.

## Limitations and troubleshooting

The command reads its positional arguments and endpoint data while importing `config.py`. It does not have an `argparse` help/validation layer: running `interp_main.py --help` or omitting arguments is not a supported discovery command. Use the four positional arguments shown above.

The printed `dQ` value comes from `get_Q`: it is a mass-weighted atomic displacement diagnostic between the endpoints, not a path energy or activation barrier.

| Symptom or concern | What to check |
| --- | --- |
| Input file cannot be found | Relative input paths start at the parent of the source checkout, not necessarily the current working directory. Follow the quick-start layout or use absolute paths. |
| `KeyError` for an element or incorrect species labels | The default configuration expects `Pb I C N H`. Do not pass another composition or order without adapting and checking the assembly code. |
| Too many or too few atoms in an extracted molecule | Inspect the cutoff, periodic images, and each `center_neighbor` file. The neighbor search is based on distance to a reference atom, not chemical bonding. |
| A rerun appears to use old molecular positions | Existing `center_neighbor/*.dat` files are reused. Start a fresh run directory with fresh endpoint copies. |
| Unexpected rotation, `NaN`, or division by zero | Check the selected reference atoms, nearly collinear vectors, spherical-coordinate poles, and angle branch crossings. This angle-based implementation has singular cases. |
| Last image differs from the relaxed final molecule | The initial internal geometry is retained. Internal-coordinate interpolation is not implemented. |
| Molecular distances change for different endpoint cells | The interpolated lattice adds affine deformation to the molecular coordinates. Rigid-molecule preservation is conditional on a fixed cell. |
| Molecules overlap or approach the framework too closely | Geometry interpolation has no energy or force evaluation and no collision-avoidance optimization. Inspect the path before NEB relaxation. |

Periodic wrapping is component-wise in fractional coordinates, including in the molecular neighbor search. It is not a general closest-image search for strongly skewed cells. Large displacements, skewed cells, symmetry-related atom permutations, and large rotations require particular care.

## Verification

The 2026-10-01 review followed the complete command pipeline, including periodic adjustment, molecular extraction, angle reconstruction, framework/lattice interpolation, POSCAR assembly, and the final coordinate shift. A new synthetic 12-atom fixed-cell rotation produced six finite, readable structures with the expected cell and species counts. The largest change in an intramolecular pair distance relative to the first generated image was approximately **2.1e-4 angstrom**, consistent with the implementation's finite-precision transformations for that fixture. This is a workflow check, not a general error bound or an NEB convergence result.

After running your own example, this read-only check validates the numbered outputs from the run directory:

```python
import numpy as np
from pymatgen.core import Structure

n_images = 6  # Total structures, including endpoints
images = [Structure.from_file(f"0{i}/POSCAR") for i in range(n_images)]
reference_species = [str(site.specie) for site in images[0]]
for image in images:
    assert [str(site.specie) for site in image] == reference_species
    assert np.isfinite(image.cart_coords).all()
    assert np.isfinite(image.lattice.matrix).all()
    assert image.volume > 0
print(f"Read {len(images)} structures with {len(images[0])} atoms each")
```

Passing these checks does not detect all short contacts, wrong molecular membership, or undesired rotational branches. Complete the physical inspection in the quick-start workflow before launching an electronic-structure calculation.

## Code map

| File | Responsibility |
| --- | --- |
| [`interp_main.py`](interp_main.py) | Command entry point. |
| [`config.py`](config.py) | Input arguments, working path, species, masses, and molecule metadata. |
| [`src/function.py`](src/function.py) | Runs the processing stages in order. |
| [`utils/mv_boundry.py`](utils/mv_boundry.py) | Rewrites the initial structure to adjust periodic displacements. |
| [`src/center.py`](src/center.py), [`src/pboundry.py`](src/pboundry.py) | Molecular grouping, periodic neighbor coordinates, and reference centers. |
| [`src/operate/write_molecule.py`](src/operate/write_molecule.py) | Writes extracted molecular data. |
| [`src/molecule_code.py`](src/molecule_code.py), [`src/coordinates_convert.py`](src/coordinates_convert.py) | Encodes molecular frames and reconstructs atomic coordinates. |
| [`src/interp_process.py`](src/interp_process.py) | Interpolates molecular centers and angles using the initial internal coordinates. |
| [`src/operate/linear_interpolate.py`](src/operate/linear_interpolate.py) | Fractional-coordinate and lattice interpolation. |
| [`src/operate/combine_poscar.py`](src/operate/combine_poscar.py) | Combines framework and molecular coordinates and writes numbered directories. |
| [`utils/reset_lattice.py`](utils/reset_lattice.py) | Applies the final coordinate translation to `interp_result/` only. |

## Contact

Xuan-Yan Chen: [xchen565@connect.hkust-gz.edu.cn](mailto:xchen565@connect.hkust-gz.edu.cn).

Report reproducible problems through [GitHub Issues](https://github.com/ShaneLogic/nebmake-mol/issues), including the command, dependency versions, and a shareable minimal endpoint pair.

## License

This project is distributed under the [GNU General Public License v3.0](LICENSE).
