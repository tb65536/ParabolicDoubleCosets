package parabolicdoublecosets;

import java.math.BigInteger;

public class ParabolicDoubleCosets {

    public static void main(String[] args) {
        long time = System.currentTimeMillis();
        int N = 100;
        System.out.println("Starting computation up to N=" + N);
        BigInteger[] fact = new BigInteger[N + 1];//factorials
        fact[0] = BigInteger.ONE;
        for (int i = 1; i <= N; i++) {
            fact[i] = fact[i - 1].multiply(BigInteger.valueOf(i));
        }
        BigInteger[][] binom = new BigInteger[N + 1][N + 1];//binomial coefficients
        for (int n = 0; n <= N; n++) {
            binom[n][0] = BigInteger.ONE;
            binom[n][n] = BigInteger.ONE;
            for (int k = 1; k < n; k++) {
                binom[n][k] = binom[n - 1][k].add(binom[n - 1][k - 1]);
            }
        }
        BigInteger[][] stirling1 = new BigInteger[N + 1][N + 1];//unsigned stirling numbers of the first kind
        for (int n = 0; n <= N; n++) {
            stirling1[n][0] = BigInteger.ZERO;
            stirling1[n][n] = BigInteger.ONE;
            for (int k = 1; k < n; k++) {
                stirling1[n][k] = stirling1[n - 1][k].multiply(BigInteger.valueOf(n - 1)).add(stirling1[n - 1][k - 1]);
            }
        }
        BigInteger[][] stirling2 = new BigInteger[N + 1][N + 1];//stirling numbers of the second kind
        for (int n = 0; n <= N; n++) {
            stirling2[n][0] = BigInteger.ZERO;
            stirling2[n][n] = BigInteger.ONE;
            for (int k = 1; k < n; k++) {
                stirling2[n][k] = stirling2[n - 1][k].multiply(BigInteger.valueOf(k)).add(stirling2[n - 1][k - 1]);
            }
        }
        BigInteger[][] f = new BigInteger[N + 1][N + 1];//generalized fubini numbers
        for (int n = 0; n <= N; n++) {
            for (int k = 0; k <= n; k++) {
                f[n][k] = BigInteger.ZERO;
                for (int j = k; j <= n; j++) {
                    f[n][k] = f[n][k].add(fact[j].multiply(stirling2[n - k][j - k]));
                }
                f[n][k] = f[n][k].divide(fact[k]);//division is expensive so we do this division ahead of time
            }
        }
        BigInteger[][] g = new BigInteger[N + 1][N + 1];//g sequence defined in the paper
        for (int n = 0; n <= N; n++) {
            for (int c = 0; c <= n; c++) {
                g[n][c] = BigInteger.ZERO;
                for (int j = 0; j <= c; j++) {
                    g[n][c] = g[n][c].add(f[n][j].multiply(f[n][c - j]));
                }
            }
        }
        BigInteger[][] T = new BigInteger[N + 1][N + 1];//T sequence defined in the paper
        for (int n = 0; n <= N; n++) {
            T[n][0] = BigInteger.ZERO;
            T[n][n] = BigInteger.ONE;
            for (int k = 1; k <= n - 1; k++) {
                T[n][k] = T[n - 1][k].multiply(BigInteger.valueOf(k * k)).add(T[n - 1][k - 1]);
            }
        }
        BigInteger[][] h = new BigInteger[N + 1][N + 1];//h sequence defined in the paper
        for (int t = 0; t <= N; t += 2) {
            for (int c = 0; 2 * c <= t; c++) {
                h[t][c] = T[t / 2][c].multiply(fact[2 * c]).multiply(BigInteger.valueOf(-1).pow(c));
            }
        }
        BigInteger[] q = new BigInteger[N + 1];//q sequence defined in the paper
        for (int n = 0; n <= N; n++) {
            q[n] = BigInteger.ZERO;
            for (int c = 0; 2 * c <= n; c++) {
                for (int t = 2 * c; t <= n; t += 2) {
                    q[n] = q[n].add(binom[n][t].multiply(h[t][c]).multiply(g[n - t + c][c]));
                }
            }
        }
        BigInteger[] p = new BigInteger[N + 1];//number of parabolic double cosets in S_n
        for (int n = 0; n <= N; n++) {
            p[n] = BigInteger.ZERO;
            for (int k = 0; k <= n; k++) {
                p[n] = p[n].add(stirling1[n][k].multiply(q[k]));
            }
            p[n] = p[n].divide(fact[n]);
        }
        System.out.println("Finished computation in " + (System.currentTimeMillis() - time) + "ms");
        for (int n = 0; n <= N; n++) {
            System.out.println(n + " " + p[n]);
        }
    }
}
