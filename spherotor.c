#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_ATOMS 2000
#define MAX_LABEL 8
#define MAX_PATH 1024
#define MAX_FRAMES 100000
#define MAX_COORDS 1000000

typedef struct {
    char el[MAX_LABEL];
    double x, y, z;
} Atom;

typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    double w, x, y, z;
} Quat;

/* ============================================================================
 * Comprehensive periodic table atomic masses
 * ============================================================================ */
static double element_mass(const char *el) {
    if (strcmp(el, "H") == 0) return 1.00794;
    if (strcmp(el, "He") == 0) return 4.00260;
    if (strcmp(el, "Li") == 0) return 6.94100;
    if (strcmp(el, "Be") == 0) return 9.01218;
    if (strcmp(el, "B") == 0) return 10.81100;
    if (strcmp(el, "C") == 0) return 12.01070;
    if (strcmp(el, "N") == 0) return 14.00670;
    if (strcmp(el, "O") == 0) return 15.99900;
    if (strcmp(el, "F") == 0) return 18.99840;
    if (strcmp(el, "Ne") == 0) return 20.17970;
    if (strcmp(el, "Na") == 0) return 22.98977;
    if (strcmp(el, "Mg") == 0) return 24.30500;
    if (strcmp(el, "Al") == 0) return 26.98154;
    if (strcmp(el, "Si") == 0) return 28.08550;
    if (strcmp(el, "P") == 0) return 30.97376;
    if (strcmp(el, "S") == 0) return 32.06500;
    if (strcmp(el, "Cl") == 0) return 35.45300;
    if (strcmp(el, "Ar") == 0) return 39.94800;
    if (strcmp(el, "K") == 0) return 39.09830;
    if (strcmp(el, "Ca") == 0) return 40.07800;
    if (strcmp(el, "Sc") == 0) return 44.95591;
    if (strcmp(el, "Ti") == 0) return 47.86700;
    if (strcmp(el, "V") == 0) return 50.94150;
    if (strcmp(el, "Cr") == 0) return 51.99610;
    if (strcmp(el, "Mn") == 0) return 54.93805;
    if (strcmp(el, "Fe") == 0) return 55.84500;
    if (strcmp(el, "Co") == 0) return 58.93320;
    if (strcmp(el, "Ni") == 0) return 58.69340;
    if (strcmp(el, "Cu") == 0) return 63.54600;
    if (strcmp(el, "Zn") == 0) return 65.39000;
    if (strcmp(el, "Ga") == 0) return 69.72300;
    if (strcmp(el, "Ge") == 0) return 72.64000;
    if (strcmp(el, "As") == 0) return 74.92160;
    if (strcmp(el, "Se") == 0) return 78.96000;
    if (strcmp(el, "Br") == 0) return 79.90400;
    if (strcmp(el, "Kr") == 0) return 83.79800;
    if (strcmp(el, "Rb") == 0) return 85.46780;
    if (strcmp(el, "Sr") == 0) return 87.62000;
    if (strcmp(el, "Y") == 0) return 88.90585;
    if (strcmp(el, "Zr") == 0) return 91.22400;
    if (strcmp(el, "Nb") == 0) return 92.90638;
    if (strcmp(el, "Mo") == 0) return 95.94000;
    if (strcmp(el, "Tc") == 0) return 98.00000;
    if (strcmp(el, "Ru") == 0) return 101.07000;
    if (strcmp(el, "Rh") == 0) return 102.90550;
    if (strcmp(el, "Pd") == 0) return 106.42000;
    if (strcmp(el, "Ag") == 0) return 107.86800;
    if (strcmp(el, "Cd") == 0) return 112.41100;
    if (strcmp(el, "In") == 0) return 114.81800;
    if (strcmp(el, "Sn") == 0) return 118.71000;
    if (strcmp(el, "Sb") == 0) return 121.75700;
    if (strcmp(el, "Te") == 0) return 127.60000;
    if (strcmp(el, "I") == 0) return 126.90447;
    if (strcmp(el, "Xe") == 0) return 131.29300;
    if (strcmp(el, "Cs") == 0) return 132.90545;
    if (strcmp(el, "Ba") == 0) return 137.32700;
    if (strcmp(el, "La") == 0) return 138.90547;
    if (strcmp(el, "Ce") == 0) return 140.11600;
    if (strcmp(el, "Pr") == 0) return 140.90765;
    if (strcmp(el, "Nd") == 0) return 144.24200;
    if (strcmp(el, "Pm") == 0) return 145.00000;
    if (strcmp(el, "Sm") == 0) return 150.36000;
    if (strcmp(el, "Eu") == 0) return 151.96400;
    if (strcmp(el, "Gd") == 0) return 157.25000;
    if (strcmp(el, "Tb") == 0) return 158.92534;
    if (strcmp(el, "Dy") == 0) return 162.50000;
    if (strcmp(el, "Ho") == 0) return 164.93032;
    if (strcmp(el, "Er") == 0) return 167.25900;
    if (strcmp(el, "Tm") == 0) return 168.93421;
    if (strcmp(el, "Yb") == 0) return 173.04000;
    if (strcmp(el, "Lu") == 0) return 174.96700;
    if (strcmp(el, "Hf") == 0) return 178.49200;
    if (strcmp(el, "Ta") == 0) return 180.94788;
    if (strcmp(el, "W") == 0) return 183.84000;
    if (strcmp(el, "Re") == 0) return 186.20700;
    if (strcmp(el, "Os") == 0) return 190.23000;
    if (strcmp(el, "Ir") == 0) return 192.21700;
    if (strcmp(el, "Pt") == 0) return 195.08400;
    if (strcmp(el, "Au") == 0) return 196.96657;
    if (strcmp(el, "Hg") == 0) return 200.59200;
    if (strcmp(el, "Tl") == 0) return 204.38300;
    if (strcmp(el, "Pb") == 0) return 207.20000;
    if (strcmp(el, "Bi") == 0) return 208.98038;
    if (strcmp(el, "Po") == 0) return 209.00000;
    if (strcmp(el, "At") == 0) return 210.00000;
    if (strcmp(el, "Rn") == 0) return 222.00000;
    if (strcmp(el, "Fr") == 0) return 223.00000;
    if (strcmp(el, "Ra") == 0) return 226.00000;
    if (strcmp(el, "Ac") == 0) return 227.00000;
    if (strcmp(el, "Th") == 0) return 232.03806;
    if (strcmp(el, "Pa") == 0) return 231.03588;
    if (strcmp(el, "U") == 0) return 238.02891;
    if (strcmp(el, "Np") == 0) return 237.00000;
    if (strcmp(el, "Pu") == 0) return 244.00000;
    if (strcmp(el, "Am") == 0) return 243.00000;
    if (strcmp(el, "Cm") == 0) return 247.00000;
    if (strcmp(el, "Bk") == 0) return 247.00000;
    if (strcmp(el, "Cf") == 0) return 251.00000;
    if (strcmp(el, "Es") == 0) return 252.00000;
    if (strcmp(el, "Fm") == 0) return 257.00000;
    return 12.0;
}

/* ============================================================================
 * Vec3 utilities
 * ============================================================================ */
static Vec3 v_sub(Vec3 a, Vec3 b) { Vec3 c = {a.x - b.x, a.y - b.y, a.z - b.z}; return c; }
static Vec3 v_scale(Vec3 a, double s) { Vec3 c = {a.x * s, a.y * s, a.z * s}; return c; }

/* ============================================================================
 * Quaternion operations
 * ============================================================================ */

/* Rotate vector v by unit quaternion q (optimized qvq* form) */
static Vec3 qrot(Quat q, Vec3 v) {
    double tx = 2.0 * (q.y * v.z - q.z * v.y);
    double ty = 2.0 * (q.z * v.x - q.x * v.z);
    double tz = 2.0 * (q.x * v.y - q.y * v.x);
    Vec3 r;
    r.x = v.x + q.w * tx + (q.y * tz - q.z * ty);
    r.y = v.y + q.w * ty + (q.z * tx - q.x * tz);
    r.z = v.z + q.w * tz + (q.x * ty - q.y * tx);
    return r;
}

/* ============================================================================
 * Uniform SO(3) sampling via R3 quasi-random sequence + Shoemake parameterization
 *
 * Shoemake (1992) showed that uniform rotations can be parameterized by three
 * independent uniform variates (s, sigma, tau) in [0,1) x [0,2pi) x [0,2pi):
 *
 *   q = ( sqrt(1-s)*sin(sigma), sqrt(1-s)*cos(sigma),
 *         sqrt(s)*sin(tau),     sqrt(s)*cos(tau)     )
 *
 * For low-discrepancy (quasi-random) coverage, we replace (s, sigma, tau) with
 * the R3 sequence -- the 3D generalization of the golden-ratio sequence.
 *
 * The R3 base constant is the real root of x^4 = x + 1:
 *   phi3 ~ 1.2207440846057596
 *   alpha1 = 1/phi3,  alpha2 = 1/phi3^2,  alpha3 = 1/phi3^3
 *
 * For sample i:
 *   s_i     = frac(0.5 + i * alpha1)
 *   sigma_i = 2*pi * frac(0.5 + i * alpha2)
 *   tau_i   = 2*pi * frac(0.5 + i * alpha3)
 *
 * This gives excellent uniformity on SO(3) for any N, with no need to
 * decompose N into factors, and every atom's sphere gets equal coverage.
 * ============================================================================ */

/* R3 quasi-random constants (real root of x^4 - x - 1 = 0) */
#define R3_PHI  1.2207440846057596
#define R3_A1   (1.0 / R3_PHI)                    /* ~0.8191725134 */
#define R3_A2   (1.0 / (R3_PHI * R3_PHI))         /* ~0.6710436067 */
#define R3_A3   (1.0 / (R3_PHI * R3_PHI * R3_PHI)) /* ~0.5497004779 */

static Quat uniform_so3_sample(int i) {
    double s     = fmod(0.5 + (double)i * R3_A1, 1.0);
    double sigma = 2.0 * M_PI * fmod(0.5 + (double)i * R3_A2, 1.0);
    double tau   = 2.0 * M_PI * fmod(0.5 + (double)i * R3_A3, 1.0);

    double r1 = sqrt(1.0 - s);
    double r2 = sqrt(s);

    Quat q;
    q.w = r1 * sin(sigma);
    q.x = r1 * cos(sigma);
    q.y = r2 * sin(tau);
    q.z = r2 * cos(tau);
    return q;
}

/* ============================================================================
 * XYZ I/O
 * ============================================================================ */
static int read_xyz(const char *path, Atom *atoms, int *n_atoms, char *title, size_t title_sz) {
    FILE *fp = fopen(path, "r");
    if (!fp) return 0;

    int n = 0;
    if (fscanf(fp, "%d", &n) != 1 || n < 1 || n > MAX_ATOMS) {
        fclose(fp);
        return 0;
    }

    int ch;
    while ((ch = fgetc(fp)) != '\n' && ch != EOF) {}

    if (!fgets(title, (int)title_sz, fp)) {
        fclose(fp);
        return 0;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%7s %lf %lf %lf", atoms[i].el, &atoms[i].x, &atoms[i].y, &atoms[i].z) != 4) {
            fclose(fp);
            return 0;
        }
    }

    fclose(fp);
    *n_atoms = n;
    return 1;
}

static void write_xyz_frame_single(FILE *fp, const Atom *atoms, int n_atoms, int frame_num) {
    fprintf(fp, "%d\n", n_atoms);
    fprintf(fp, "Frame %d\n", frame_num);
    for (int i = 0; i < n_atoms; i++) {
        fprintf(fp, "%-2s % .9f % .9f % .9f\n", atoms[i].el, atoms[i].x, atoms[i].y, atoms[i].z);
    }
}

/* ============================================================================
 * Filesystem and coordinate helpers
 * ============================================================================ */
static int ensure_dir(const char *dir) {
    struct stat st;
    if (stat(dir, &st) == 0) {
        return S_ISDIR(st.st_mode) ? 1 : 0;
    }
    if (mkdir(dir, 0755) == 0) return 1;
    return (errno == EEXIST);
}

static int read_coords_txt(const char *path, Vec3 **coords_out, int *n_coords_out) {
    FILE *fp = fopen(path, "r");
    if (!fp) return 0;

    int cap = 1024;
    int n = 0;
    Vec3 *coords = (Vec3 *)malloc((size_t)cap * sizeof(Vec3));
    if (!coords) {
        fclose(fp);
        return 0;
    }

    char line[512];
    while (fgets(line, (int)sizeof(line), fp)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '#') continue;

        double x, y, z;
        char lbl[64];
        int m = sscanf(p, "%lf %lf %lf", &x, &y, &z);
        if (m != 3) {
            m = sscanf(p, "%63s %lf %lf %lf", lbl, &x, &y, &z);
        }
        if (m != 3 && m != 4) continue;

        if (n >= MAX_COORDS) {
            free(coords);
            fclose(fp);
            return 0;
        }

        if (n == cap) {
            int new_cap = cap * 2;
            Vec3 *tmp = (Vec3 *)realloc(coords, (size_t)new_cap * sizeof(Vec3));
            if (!tmp) {
                free(coords);
                fclose(fp);
                return 0;
            }
            coords = tmp;
            cap = new_cap;
        }

        coords[n].x = x;
        coords[n].y = y;
        coords[n].z = z;
        n++;
    }

    fclose(fp);
    if (n < 1) {
        free(coords);
        return 0;
    }

    *coords_out = coords;
    *n_coords_out = n;
    return 1;
}

/* ============================================================================
 * Core rotation loop -- uniform SO(3) via R3 quasi-random quaternions
 * ============================================================================ */
static int run(const char *in_path, const char *out_dir, const char *prefix,
               int n_frames, const char *coords_path) {

    Atom base_atoms[MAX_ATOMS];
    int n_atoms = 0;
    char in_title[256];

    if (!read_xyz(in_path, base_atoms, &n_atoms, in_title, sizeof(in_title))) {
        fprintf(stderr, "Error: failed reading XYZ: %s\n", in_path);
        return 1;
    }

    /* Compute center of mass */
    Vec3 com = {0.0, 0.0, 0.0};
    double msum = 0.0;
    for (int i = 0; i < n_atoms; i++) {
        double m = element_mass(base_atoms[i].el);
        com.x += m * base_atoms[i].x;
        com.y += m * base_atoms[i].y;
        com.z += m * base_atoms[i].z;
        msum += m;
    }
    if (msum <= 0.0) {
        fprintf(stderr, "Error: invalid mass sum.\n");
        return 1;
    }
    com = v_scale(com, 1.0 / msum);

    /* Allocate working arrays */
    Vec3 *com_frame = (Vec3 *)malloc((size_t)n_atoms * sizeof(Vec3));
    Atom *out_atoms = (Atom *)malloc((size_t)n_atoms * sizeof(Atom));
    if (!com_frame || !out_atoms) {
        fprintf(stderr, "Error: memory allocation failed.\n");
        free(com_frame);
        free(out_atoms);
        return 1;
    }

    /* Shift atoms to COM frame */
    for (int i = 0; i < n_atoms; i++) {
        Vec3 r = {base_atoms[i].x, base_atoms[i].y, base_atoms[i].z};
        com_frame[i] = v_sub(r, com);
        strcpy(out_atoms[i].el, base_atoms[i].el);
    }

    /* Load translation centers (optional) */
    Vec3 *centers = NULL;
    int n_centers = 0;
    int has_coords = 0;
    if (coords_path) {
        if (!read_coords_txt(coords_path, &centers, &n_centers)) {
            fprintf(stderr, "Error: failed reading coordinate file: %s\n", coords_path);
            free(com_frame);
            free(out_atoms);
            return 1;
        }
        has_coords = 1;
    } else {
        centers = (Vec3 *)malloc(sizeof(Vec3));
        if (!centers) {
            fprintf(stderr, "Error: memory allocation failed.\n");
            free(com_frame);
            free(out_atoms);
            return 1;
        }
        centers[0] = com;
        n_centers = 1;
    }

    /* Open trajectory output */
    char traj_path[MAX_PATH];
    snprintf(traj_path, sizeof(traj_path), "%s/%s.xyz", out_dir, prefix);

    FILE *traj_fp = fopen(traj_path, "w");
    if (!traj_fp) {
        fprintf(stderr, "Error: cannot open trajectory file: %s\n", traj_path);
        free(com_frame);
        free(out_atoms);
        free(centers);
        return 1;
    }

    /* For each frame, generate a single quaternion from the R3 quasi-random
     * sequence on SO(3). Apply it to all atoms. Every atom gets uniformly
     * distributed sphere coverage with no axis dependence. */
    int written = 0;
    int sample_idx = 0;

    for (int cidx = 0; cidx < n_centers; cidx++) {
        Vec3 center = centers[cidx];

        for (int f = 0; f < n_frames; f++) {
            Quat q = uniform_so3_sample(sample_idx++);

            for (int i = 0; i < n_atoms; i++) {
                Vec3 rotated = qrot(q, com_frame[i]);
                out_atoms[i].x = rotated.x + center.x;
                out_atoms[i].y = rotated.y + center.y;
                out_atoms[i].z = rotated.z + center.z;
            }

            write_xyz_frame_single(traj_fp, out_atoms, n_atoms, written + 1);
            written++;
        }
    }

    fclose(traj_fp);
    if (has_coords) {
        printf("Done. Wrote %d frames (%d centers x %d uniform SO(3) samples) to: %s\n",
               written, n_centers, n_frames, traj_path);
    } else {
        printf("Done. Wrote %d frames (%d uniform SO(3) samples) to: %s\n",
               written, n_frames, traj_path);
    }

    free(com_frame);
    free(out_atoms);
    free(centers);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 5 || argc > 6) {
        fprintf(stderr,
                "Usage:\n"
                "  %s input.xyz out_dir out_prefix n_frames [coords.txt]\n\n"
                "Description:\n"
                "  Generates uniformly sampled molecular rotations via R3\n"
                "  quasi-random quaternions (Shoemake parameterization of SO(3)).\n"
                "  Every atom traces a uniformly covered sphere -- no axis\n"
                "  selection needed, no product-grid artifacts.\n\n"
                "  If coords.txt is provided, the rotated molecule is placed\n"
                "  at each coordinate, yielding centers x n_frames total frames.\n\n"
                "  Full periodic table supported (H through Fm).\n",
                argv[0]);
        return 1;
    }

    int n_frames = atoi(argv[4]);
    const char *coords_path = (argc == 6) ? argv[5] : NULL;

    if (n_frames < 1 || n_frames > MAX_FRAMES) {
        fprintf(stderr, "Error: n_frames must be between 1 and %d.\n", MAX_FRAMES);
        return 1;
    }
    if (!ensure_dir(argv[2])) {
        fprintf(stderr, "Error: cannot create/access output directory: %s\n", argv[2]);
        return 1;
    }

    return run(argv[1], argv[2], argv[3], n_frames, coords_path);
}



// ./spherotor input.xyz out_dir out_prefix n_frames [coords.txt]
// Example:
//   ./spherotor molecule.xyz output mol 1000
//   ./spherotor molecule.xyz output mol 1000 coords.txt