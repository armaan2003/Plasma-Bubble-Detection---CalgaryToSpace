// Computes the S4 amplitude scintillation index from C/N0 for every satellite signal in 60 s windows (Kintner 2007).
#include <stdio.h>
#include <math.h>

#define MAX 400000

int n;
double t[MAX], cn0[MAX], st[MAX], sc[MAX];
int sys[MAX], prn[MAX], sig[MAX];
char done[8][256][32];

// Loads every phase locked row from the parser CSV.
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
        int ph;
        if (sscanf(line, "%lf,%*d,%d,%d,%d,%lf,%*f,%*f,%*f,%*f,%d",
                   &t[n], &sys[n], &prn[n], &sig[n], &cn0[n], &ph) == 6 && ph == 1 && prn[n] < 256)
            n++;
    }
    fclose(f);
}

double s4(double *c, int m)
{
    double sum = 0, sum_sq = 0;
    for (int k = 0; k < m; k++) {
        double I = pow(10, c[k] / 10);
        sum += I;
        sum_sq += I * I;
    }
    double mean = sum / m;
    double var = sum_sq / m - mean * mean;
    // Rounding can make the variance slightly negative.
    if (var < 0) var = 0;
    return sqrt(var) / mean;
}

int main(int argc, char **argv)
{
    char *file = "lab_5hz.csv";
    if (argc > 1) file = argv[1];
    load(file);
    printf("%d rows from %s\n", n, file);
    printf("sys prn sig   start_s    samples   S4\n");

    for (int i = 0; i < n; i++) {
        if (done[sys[i]][prn[i]][sig[i]]) continue;
        done[sys[i]][prn[i]][sig[i]] = 1;

        int m = 0;
        for (int k = i; k < n; k++) {
            if (sys[k] == sys[i] && prn[k] == prn[i] && sig[k] == sig[i]) {
                st[m] = t[k];
                sc[m] = cn0[k];
                m++;
            }
        }

        // The last unfinished window is dropped.
        int start = 0;
        for (int k = 0; k < m; k++) {
            if (st[k] - st[start] >= 60) {
                double s = s4(&sc[start], k - start);
                char *level = "";
                if (s > 0.1) level = "weak";
                if (s > 0.5) level = "STRONG";
                printf("%3d %3d %3d  %.1f  %6d    %.3f %s\n", sys[i], prn[i], sig[i], st[start], k - start, s, level);
                start = k;
            }
        }
    }
    return 0;
}
