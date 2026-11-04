// Computes the sigma phi phase scintillation index from L1 carrier phase for every satellite in 60 s windows (Kintner 2007).
#include <stdio.h>
#include <math.h>

#define MAX 400000
#define PI 3.14159265358979
// RANGECMP ADR rolls over every 8388608 cycles.
#define WRAP 8388608.0

int n;
double t[MAX], cn0[MAX], adr[MAX], lock[MAX];
int sys[MAX], prn[MAX];
int nsat, list_sys[256], list_prn[256];
double rt[MAX], ra[MAX], rl[MAX];
double pt[MAX], pa[MAX], pl[MAX];
double xt[MAX], xa[MAX];
char slip[MAX];

// Loads GPS L1 C/A, Galileo E1 and QZSS L1 C/A rows with phase lock and known parity.
void load(char *file)
{
    FILE *f = fopen(file, "r");
    if (!f) {
        printf("cannot open %s\n", file);
        return;
    }
    char line[512];
    fgets(line, sizeof line, f);
    while (n < MAX && fgets(line, sizeof line, f)) {
        int g, ph, par;
        if (sscanf(line, "%lf,%*d,%d,%d,%d,%lf,%lf,%*f,%*f,%lf,%d,%d",
                   &t[n], &sys[n], &prn[n], &g, &cn0[n], &adr[n], &lock[n], &ph, &par) != 9) continue;
        int l1 = (sys[n] == 0 && g == 0) || (sys[n] == 3 && g == 2) || (sys[n] == 5 && g == 0);
        if (l1 && ph == 1 && par == 1) n++;
    }
    fclose(f);
}

// Copies one satellite into the given arrays and undoes the RANGECMP rollover.
int pick(int s, int p, double *to, double *ao, double *lo)
{
    int m = 0;
    for (int k = 0; k < n; k++) {
        if (sys[k] != s || prn[k] != p) continue;
        to[m] = t[k];
        ao[m] = adr[k];
        lo[m] = lock[k];
        if (m > 0 && ao[m] - ao[m - 1] > WRAP / 2) ao[m] -= WRAP;
        if (m > 0 && ao[m] - ao[m - 1] < -WRAP / 2) ao[m] += WRAP;
        m++;
    }
    return m;
}

// Detrends with a straight line fitted over 10 s because samples are unevenly spaced and the phase ramps fast.
double sigma_phi(double *tm, double *a, int cnt)
{
    int count = 0;
    double sum = 0, sum_sq = 0;
    for (int k = 0; k < cnt; k++) {
        if (tm[k] - tm[0] < 5 || tm[cnt - 1] - tm[k] < 5) continue;

        double mt = 0, ma = 0;
        int c = 0;
        for (int j = 0; j < cnt; j++) {
            if (fabs(tm[j] - tm[k]) <= 5) {
                mt += tm[j];
                ma += a[j];
                c++;
            }
        }
        mt = mt / c;
        ma = ma / c;

        double top = 0, bottom = 0;
        for (int j = 0; j < cnt; j++) {
            if (fabs(tm[j] - tm[k]) <= 5) {
                top += (tm[j] - mt) * (a[j] - ma);
                bottom += (tm[j] - mt) * (tm[j] - mt);
            }
        }
        double trend = ma + (top / bottom) * (tm[k] - mt);

        double phi = 2 * PI * (a[k] - trend);
        sum += phi;
        sum_sq += phi * phi;
        count++;
    }
    double mean = sum / count;
    double var = sum_sq / count - mean * mean;
    // Rounding can make the variance slightly negative.
    if (var < 0) var = 0;
    return sqrt(var);
}

int main(int argc, char **argv)
{
    char *file = "lab_5hz.csv";
    if (argc > 1) file = argv[1];
    load(file);
    printf("%d L1 rows from %s\n", n, file);

    for (int i = 0; i < n; i++) {
        int found = 0;
        for (int j = 0; j < nsat; j++)
            if (list_sys[j] == sys[i] && list_prn[j] == prn[i]) found = 1;
        if (!found && nsat < 256) {
            list_sys[nsat] = sys[i];
            list_prn[nsat] = prn[i];
            nsat++;
        }
    }

    // The reference is the satellite with the highest mean C/N0.
    double best = 0;
    int rsys = 0, rprn = 0;
    for (int j = 0; j < nsat; j++) {
        double sum = 0;
        int c = 0;
        for (int k = 0; k < n; k++) {
            if (sys[k] == list_sys[j] && prn[k] == list_prn[j]) {
                sum += cn0[k];
                c++;
            }
        }
        if (sum / c > best) {
            best = sum / c;
            rsys = list_sys[j];
            rprn = list_prn[j];
        }
    }
    int nr = pick(rsys, rprn, rt, ra, rl);
    printf("reference: sys %d prn %d, mean C/N0 %.1f dB-Hz\n", rsys, rprn, best);
    printf("sys prn   start_s    samples   sigma_phi (rad)\n");

    for (int j = 0; j < nsat; j++) {
        if (list_sys[j] == rsys && list_prn[j] == rprn) continue;

        // Subtracting the reference at matching times removes the receiver clock.
        int np = pick(list_sys[j], list_prn[j], pt, pa, pl);
        int m = 0, r = 0;
        for (int k = 0; k < np; k++) {
            while (r < nr && rt[r] < pt[k]) r++;
            if (r == nr) break;
            if (rt[r] != pt[k]) continue;
            xt[m] = pt[k];
            xa[m] = pa[k] - ra[r];
            // Lock time going down means a cycle slip.
            slip[m] = m > 0 && (pl[k] < pl[k - 1] || rl[r] < rl[r - 1]);
            m++;
        }

        // The last unfinished window is dropped.
        int start = 0;
        for (int k = 0; k < m; k++) {
            if (xt[k] - xt[start] < 60) continue;

            // Windows with a gap over 1 s or a cycle slip are skipped.
            int ok = 1;
            for (int i = start + 1; i < k; i++)
                if (xt[i] - xt[i - 1] > 1 || slip[i]) ok = 0;

            printf("%3d %3d  %.1f  %6d    ", list_sys[j], list_prn[j], xt[start], k - start);
            if (ok)
                printf("%.3f\n", sigma_phi(&xt[start], &xa[start], k - start));
            else
                printf("skipped (gap or cycle slip)\n");
            start = k;
        }
    }
    return 0;
}
