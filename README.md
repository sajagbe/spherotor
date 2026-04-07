# spherotor

A command-line tool that generates a **multi-frame XYZ trajectory** by rotating a single molecule through uniformly distributed orientations on a sphere.

## Basic Usage

```bash
./spherotor input.xyz out_dir out_prefix n_points axis_i axis_j [psi_deg] [coords.txt]
```

---

## What It Does

1. Reads a molecule from `input.xyz` and computes its **center of mass (COM)**
2. Defines a **molecular axis** as the vector from atom `axis_i` → `axis_j`
3. Generates `n_points` roughly-uniform directions on the sphere (Fibonacci method)
4. For each direction, rotates the molecule so the molecular axis aligns to it
5. Writes all frames to `out_dir/out_prefix.xyz`

---

## Modes

| Mode | Command | Frames | Behavior |
| --- | --- | --- | --- |
| **Spherical rotation only** | No optional args | `n_points` | Molecule stays at original COM, rotated through all directions |
| **With translation file** | `coords.txt` provided | `n_centers × n_points` | Full rotation set repeated at each coordinate center |
| **With twist** | `psi_deg` provided | same as above | Adds a rotation about the target axis after alignment |

> Frame ordering when using a coords file: all `n_points` rotations at center 1, then center 2, etc.
> 

---

## Optional Argument Parsing

- **1 optional arg after `axis_j`** → number = `psi_deg`; non-number = `coords.txt`
- **2 optional args** → first is `psi_deg`, second is `coords.txt`

## To Compile:

```bash
cc -O3 -std=c11 -Wall -Wextra -pedantic -o spherotor spherotor.c -lm 
```

e.g

To translate to origin alone (no COM rotation, no psi rotation)

```bash
./spherotor LF.xyz . LF-spherotor 1 1 2 coord.txt
```

Here `coord.txt` contains just `0.000 0.000 0.000` .

Next I copy this to make one where the atoms are placed by placeholder X. In vsg, atoms named X  are treated with vdw radius of 1Å.

Then:

```bash
vsg -s 7 -d 0.04 LF-X.xyz -t 
vsg -s 9 -d 0.04 LF-X.xyz -t 
```

0.04 density because we want to use 1 COM coordinate per 25 $Å^2$, we chose this because we wanted to minimize crossover of sampled area for each copy.  7 and 9 Å away because we wanna minimize vdW interaction with LF. 

```bash
./spherotor p-cresol.xyz . p-cresol-7 50 1 2 LF-X7.txt
Done. Wrote 2200 frames (44 centers x 50 rotations) to: ./p-cresol-7.xyz

./spherotor p-cresol.xyz . p-cresol-9 50 1 2 LF-X9.txt
Done. Wrote 3300 frames (66 centers x 50 rotations) to: ./p-cresol-9.xyz

 ./spherotor Indole.xyz . Indole-9 50 1 2 LF-X9.txt
Done. Wrote 3300 frames (66 centers x 50 rotations) to: ./Indole-9.xyz

./spherotor Indole.xyz . Indole-7 50 1 2 LF-X7.txt
Done. Wrote 2200 frames (44 centers x 50 rotations) to: ./Indole-7.xyz
```