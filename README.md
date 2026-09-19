# spherotor

A command-line tool that generates a **multi-frame XYZ trajectory** by rotating a single molecule through uniformly distributed orientations on a sphere (full SO(3), not just a single axis).

## Basic Usage

```bash
./spherotor input.xyz out_dir out_prefix n_frames [coords.txt]
```

---

## What It Does

1. Reads a molecule from `input.xyz` and computes its **center of mass (COM)**.
2. Shifts the molecule into the COM frame.
3. Generates `n_frames` uniformly-distributed rotations of the *whole molecule* using a low-discrepancy R3 quasi-random sequence fed through Shoemake's quaternion parameterization of SO(3) — this covers all possible orientations evenly, with no molecular axis needed and no product-grid artifacts.
4. For each rotation, places the rotated molecule at a center (either the original COM, or each coordinate in `coords.txt` if provided).
5. Writes all frames to `out_dir/out_prefix.xyz`.

---

## Modes

| Mode | Command | Frames | Behavior |
| --- | --- | --- | --- |
| **Spherical rotation only** | No `coords.txt` | `n_frames` | Molecule stays at its original COM, rotated through `n_frames` orientations |
| **With translation file** | `coords.txt` provided | `n_centers × n_frames` | Full rotation set repeated at each coordinate center |

> Frame ordering when using a coords file: all `n_frames` rotations at center 1, then center 2, etc.

## To Compile:

```bash
cc -O3 -std=c11 -Wall -Wextra -pedantic -o spherotor spherotor.c -lm 
```

e.g

To translate to origin alone (no rotation sampling, single frame per center)

```bash
./spherotor LF.xyz . LF-spherotor 1 coord.txt
```

Here `coord.txt` contains just `0.000 0.000 0.000` .

Next I copy this to make one where the atoms are placed by placeholder X. In vsg, atoms named X are treated with vdw radius of 1Å.

Then:

```bash
vsg -s 7 -d 0.04 LF-X.xyz -t 
vsg -s 9 -d 0.04 LF-X.xyz -t 
```

0.04 density because we want to use 1 COM coordinate per 25 $Å^2$ (1/0.04 = 25), we chose this because we wanted to minimize crossover of sampled area for each copy. 7 and 9 Å away because we wanna minimize vdW interaction with LF.

```bash
./spherotor p-cresol.xyz 7 cresol7 50 LF-Prep/LF-X7.txt
Done. Wrote 2150 frames (43 centers x 50 uniform SO(3) samples) to: 7/cresol7.xyz

./spherotor p-cresol.xyz 9 cresol9 50 LF-Prep/LF-X9.txt
Done. Wrote 3250 frames (65 centers x 50 uniform SO(3) samples) to: 9/cresol9.xyz

./spherotor Indole.xyz 9 indole9 50 LF-Prep/LF-X9.txt
Done. Wrote 3250 frames (65 centers x 50 uniform SO(3) samples) to: 9/indole9.xyz

./spherotor Indole.xyz 7 indole7 50 LF-Prep/LF-X7.txt
Done. Wrote 2150 frames (43 centers x 50 uniform SO(3) samples) to: 7/indole7.xyz
```
