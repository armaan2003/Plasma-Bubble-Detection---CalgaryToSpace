// Computes ROTI for every GPS satellite with L1 and L2 using 30 s ROT and 5 min windows (Pi 1997).
#include <stdio.h>
#include <math.h>

#define MAX 400000
// GPS L2C signal, use 9 for L2P(Y).
#define L2SIG 17

int n;
double t[MAX], adr[MAX], lock[MAX];
int prn[MAX], sig[MAX];
double t1[MAX], a1[MAX], k1[MAX];
double t2[MAX], a2[MAX], k2[MAX];
double tt[MAX], tec[MAX];
char slip[MAX];
double trot[MAX], rot[MAX];

// Loads GPS L1 and L2 rows on whole 30 s epochs with phase lock and known parity.
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
        int s, ph, par;
        if (sscanf(line, "%lf,%*d,%d,%d,%d,%*f,%lf,%*f,%*f,%lf,%d,%d",
                   &t[n], &s, &prn[n], &sig[n], &adr[n], &lock[n], &ph, &par) != 8) continue;
        int epoch30 = fabs(t[n] - 30 * floor(t[n] / 30 + 0.5)) < 0.001;
        if (s == 0 && (sig[n] == 0 || sig[n] == L2SIG) && ph == 1 && par == 1 && epoch30) n++;
    }
    fclose(f);
}

int pick(int p, int g, double *to, double *ao, double *lo)
{
    int m = 0;
    for (int k = 0; k < n; k++) {
        if (prn[k] == p && sig[k] == g) {
            to[m] = t[k];
            ao[m] = adr[k];
            lo[m] = lock[k];
            m++;
        }
    }
    return m;
}

double roti(double *r, int cnt)
{
    double sum = 0, sum_sq = 0;
    for (int k = 0; k < cnt; k++) {
        sum += r[k];
        sum_sq += r[k] * r[k];
    }
    double mean = sum / cnt;
    double var = sum_sq / cnt - mean * mean;
    // Rounding can make the variance slightly negative.
    if (var < 0) var = 0;
    return sqrt(var);
}

int main(int argc, char **argv)
{
    char *file = "lab_5hz.csv";
    if (argc > 1) file = argv[1];
    load(file);
    printf("%d GPS rows on 30 s epochs from %s\n", n, file);

    double c = 299792458, f1 = 1575.42e6, f2 = 1227.60e6;
    int found = 0;
    printf("prn   start_s    ROTs   ROTI (TECU/min)\n");

    for (int p = 1; p <= 32; p++) {
        int n1 = pick(p, 0, t1, a1, k1);
        int n2 = pick(p, L2SIG, t2, a2, k2);
        if (n1 == 0 || n2 == 0) continue;
        found++;

        int m = 0, j = 0;
        for (int k = 0; k < n1; k++) {
            while (j < n2 && t2[j] < t1[k]) j++;
            if (j == n2) break;
            if (t2[j] != t1[k]) continue;
            tt[m] = t1[k];
            // NovAtel ADR has the opposite sign of carrier phase so this is L2 minus L1.
            double gf = (c / f2) * a2[j] - (c / f1) * a1[k];
            tec[m] = gf / (40.3 * (1 / (f2 * f2) - 1 / (f1 * f1))) / 1e16;
            // Lock time shorter than the step means lock was lost in between.
            slip[m] = m > 0 && (k1[k] < tt[m] - tt[m - 1] || k2[j] < tt[m] - tt[m - 1]);
            m++;
        }

        // ROT is only taken between samples exactly 30 s apart with no cycle slip.
        int nrot = 0;
        for (int k = 1; k < m; k++) {
            if (fabs(tt[k] - tt[k - 1] - 30) > 0.001 || slip[k]) continue;
            trot[nrot] = tt[k];
            rot[nrot] = (tec[k] - tec[k - 1]) / 0.5;
            nrot++;
        }

        // The last unfinished window is dropped.
        int start = 0;
        for (int k = 0; k < nrot; k++) {
            if (trot[k] - trot[start] < 300) continue;
            double r = roti(&rot[start], k - start);
            char *level = "none";
            if (r >= 0.25) level = "weak";
            if (r >= 0.5) level = "moderate";
            if (r >= 1) level = "STRONG";
            printf("%3d  %.1f  %4d    %.3f %s\n", p, trot[start], k - start, r, level);
            start = k;
        }
    }
    if (found == 0) printf("no satellite has both L1 and L2 (receiver tracks L1 only), ROTI needs L1 + L2\n");
    return 0;
}
