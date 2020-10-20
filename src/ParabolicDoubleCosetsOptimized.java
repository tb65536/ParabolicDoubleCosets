package parabolicdoublecosetsoptimized;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Arrays;

public class ParabolicDoubleCosetsOptimized {

    private static int N;
    private static BigInteger[] fact;//factorials
    
    private static final int interval = 500;//how many values of n-t+c will be run in the same batch

    public static void main(String[] args) throws IOException {
        long time = System.currentTimeMillis();
        try {
            N = Integer.parseInt(args[0]);
        } catch (Exception ex) {
            System.out.println("Error: N is required as an argument.");
            System.out.println("Usage: java -jar ParabolicDoubleCosetsOptimized.jar [N]");
            System.exit(0);
        }
        new File(statusfile).delete();
        new File(datafile).delete();
        new File(tablefile).delete();
        writeln(statusfile, "Starting computation up to N=" + N);
        writeln(tablefile, "\\[\\begin{tabular}{r|l r|l}");
        writeln(tablefile, "$n$&$p_n$&&$(p_n\\log^{2n}2)/n!$\\\\\\hline");
        fact = new BigInteger[N + 1];
        fact[0] = BigInteger.ONE;
        for (int i = 1; i <= N; i++) {
            fact[i] = fact[i - 1].multiply(BigInteger.valueOf(i));
        }
        BigInteger[] q = new BigInteger[N + 1];//q sequence
        Arrays.fill(q, BigInteger.ZERO);
        BigInteger[] stirling1 = new BigInteger[0];//unsigned stirling numbers of the first kind
        int min = 0;
        int max = interval;
        while (true) {
            if (min > N) {
                break;
            }
            if (max > N) {
                max = N;
            }
            int[] array = new int[max - min + 1];//values of n-t+c in this batch
            for (int ntc = min; ntc <= max; ntc++) {
                array[ntc - min] = ntc;
            }
            BigInteger[] qq = Arrays.stream(array).parallel().collect(() -> {//each value of n-t+c contributes to the q sequence
                BigInteger[] a = new BigInteger[N + 1];
                Arrays.fill(a, BigInteger.ZERO);
                return a;
            }, (BigInteger[] t, int value) -> {
                BigInteger[] a = q(value);
                for (int k = 0; k <= N; k++) {
                    t[k] = t[k].add(a[k]);
                }
            }, (BigInteger[] t, BigInteger[] u) -> {
                for (int k = 0; k <= N; k++) {
                    t[k] = t[k].add(u[k]);
                }
            });
            for (int k = 0; k <= N; k++) {
                q[k] = q[k].add(qq[k]);//update q squence
            }
            for (int n = min; n <= max; n++) {
                BigInteger[] newstirling1 = new BigInteger[n + 1];//update unsigned stirling numbers of the first kind
                newstirling1[0] = BigInteger.ZERO;
                newstirling1[n] = BigInteger.ONE;
                for (int k = 1; k < n; k++) {
                    newstirling1[k] = stirling1[k].multiply(BigInteger.valueOf(n - 1)).add(stirling1[k - 1]);
                }
                stirling1 = newstirling1;
                BigInteger p = BigInteger.ZERO;//number of parabolic double cosets in S_n
                for (int k = 0; k <= n; k++) {
                    p = p.add(stirling1[k].multiply(q[k]));
                }
                p = p.divide(fact[n]);
                writeln(datafile, n + " " + p);
                if (include(n)) {
                    print(n, p);
                }
            }
            writeln(statusfile, "Finished " + max + "/" + N + " after " + (System.currentTimeMillis() - time) + "ms");
            min = max + 1;
            max += interval;
        }
        writeln(tablefile, "\\end{tabular}\\]");
    }

    private static BigInteger[] q(int ntc) {
        BigInteger[] q = new BigInteger[N + 1];
        Arrays.fill(q, BigInteger.ZERO);
        BigInteger[] g = g(ntc);//g sequence for this choice of n-t+c
        BigInteger[] T = new BigInteger[0];//T sequence
        BigInteger[] binom = new BigInteger[N + 1];//binomial coefficients
        Arrays.fill(binom, BigInteger.ZERO);
        for (int t = 0; t <= N; t++) {
            BigInteger[] newB = new BigInteger[N + 1];//update binomial coefficients
            newB[t] = BigInteger.ONE;
            for (int n = t + 1; n <= N; n++) {
                newB[n] = binom[n - 1].add(newB[n - 1]);
            }
            binom = newB;
            if (t % 2 != 0) {//if t is odd then there is no contribution to q sequence
                continue;
            }
            BigInteger[] newT = new BigInteger[t / 2 + 1];//update T sequence
            newT[0] = BigInteger.ZERO;
            newT[t / 2] = BigInteger.ONE;
            for (int c = 1; c < t / 2; c++) {
                newT[c] = T[c].multiply(BigInteger.valueOf(c * c)).add(T[c - 1]);//update the T sequence
            }
            T = newT;
            int min = Math.max(0, ntc + t - N);
            int max = Math.min(t / 2, ntc);
            for (int c = min; c <= max; c++) {
                int n = ntc + t - c;
                if (c % 2 == 0) {
                    q[n] = q[n].add(binom[n].multiply(T[c]).multiply(fact[2 * c]).multiply(g[c]));
                } else {
                    q[n] = q[n].subtract(binom[n].multiply(T[c]).multiply(fact[2 * c]).multiply(g[c]));
                }
            }
        }
        return q;
    }

    private static BigInteger[] g(int n) {
        BigInteger[] f = f(n);//f sequence
        BigInteger[] g = new BigInteger[n + 1];//g sequence
        for (int c = 0; c <= n; c++) {
            g[c] = BigInteger.ZERO;
            for (int j = 0; j <= c; j++) {
                g[c] = g[c].add(f[j].multiply(f[c - j]));
            }
        }
        return g;
    }

    private static BigInteger[] f(int n) {
        BigInteger[] stirling2 = new BigInteger[0];//stirling numbers of the second kind
        BigInteger[] f = new BigInteger[n + 1];//f sequence
        for (int k = 0; k <= n; k++) {
            BigInteger[] newstirling2 = new BigInteger[k + 1];//update stirling numbers of the second kind
            newstirling2[0] = BigInteger.ZERO;
            newstirling2[k] = BigInteger.ONE;
            for (int j = 1; j < k; j++) {
                newstirling2[j] = stirling2[j].multiply(BigInteger.valueOf(j)).add(stirling2[j - 1]);
            }
            stirling2 = newstirling2;
            f[n - k] = BigInteger.ZERO;
            for (int j = 0; j <= k; j++) {
                f[n - k] = f[n - k].add(fact[n - (k - j)].multiply(stirling2[j]));
            }
            f[n - k] = f[n - k].divide(fact[n - k]);
        }
        return f;
    }

    private static boolean include(int n) {
        if (n == 0) {
            return false;
        }
        if (n <= 20) {
            return true;
        }
        if (n % 10 == 0 && n <= 100) {
            return true;
        }
        if (n % 100 == 0 && n <= 1000) {
            return true;
        }
        if (n % 500 == 0) {
            return true;
        }
        return false;
    }

    private static final int maxlen = 25;//maximum length before displaying "..."
    private static final int display = 12;//how may digits on either side of the "..."
    private static final int buffer = 12;//how many digits to use when converting from BigInteger to double

    private static void print(int n, BigInteger p) throws IOException {
        int len = p.toString().length();
        String s = p.toString();
        if (len > maxlen) {
            s = s.substring(0, display) + "..." + s.substring(len - display, len);
        }
        BigInteger b = p.multiply(BigInteger.TEN.pow(buffer)).divide(fact[n]);
        int shift = Math.max(b.toString().length() - buffer, 0);
        double d = Math.log(b.divide(BigInteger.TEN.pow(shift)).doubleValue());
        shift -= buffer;
        d += 2 * n * Math.log(Math.log(2));
        d += shift * Math.log(10);
        d = Math.exp(d);
        d = (double) Math.round(d * 1000000) / 1000000;
        writeln(tablefile, n + "&" + s + "&(" + len + " digits)&" + d + "\\\\");
    }

    private static final String statusfile = "status.txt";
    private static final String datafile = "data.txt";
    private static final String tablefile = "table.txt";

    private static void writeln(String filename, String line) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(filename, true));
        writer.write(line);
        writer.newLine();
        writer.close();
    }
}
