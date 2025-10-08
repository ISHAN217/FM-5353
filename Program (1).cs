using System;
using System.Collections.Concurrent;
using System.Threading.Tasks;

[Flags]
public enum VarianceReduction
{
    None = 0,
    Antithetic = 1 << 0,
    DeltaControlVariate = 1 << 1
}

public static class Normal
{
    // Abramowitz–Stegun 7.1.26 erf approximation
    public static double Erf(double x)
    {
        const double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741,
                     a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
        int sign = x < 0 ? -1 : 1;
        x = Math.Abs(x);
        double t = 1.0 / (1.0 + p * x);
        double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);
        return sign * y;
    }

    public static double Phi(double x) => 0.5 * (1.0 + Erf(x / Math.Sqrt(2.0)));
    public static double Pdf(double x) => Math.Exp(-0.5 * x * x) / Math.Sqrt(2.0 * Math.PI);
}

public static class BlackScholes
{
    public static double D1(double S, double K, double r, double q, double sigma, double T)
    {
        double v = sigma * Math.Sqrt(T);
        return (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / v;
    }

    public static double PriceCall(double S, double K, double r, double q, double sigma, double T)
    {
        double d1 = D1(S, K, r, q, sigma, T);
        double d2 = d1 - sigma * Math.Sqrt(T);
        return S * Math.Exp(-q * T) * Normal.Phi(d1) - K * Math.Exp(-r * T) * Normal.Phi(d2);
    }

    public static double PricePut(double S, double K, double r, double q, double sigma, double T)
    {
        double d1 = D1(S, K, r, q, sigma, T);
        double d2 = d1 - sigma * Math.Sqrt(T);
        return K * Math.Exp(-r * T) * Normal.Phi(-d2) - S * Math.Exp(-q * T) * Normal.Phi(-d1);
    }

    public static double CallDelta(double S, double K, double r, double q, double sigma, double T)
    {
        double d1 = D1(S, K, r, q, sigma, T);
        return Math.Exp(-q * T) * Normal.Phi(d1);
    }

    public static double PutDelta(double S, double K, double r, double q, double sigma, double T)
    {
        double d1 = D1(S, K, r, q, sigma, T);
        return Math.Exp(-q * T) * (Normal.Phi(d1) - 1.0);
    }
}

public interface IPathGen
{
    // Returns primary path and optional antithetic path over a grid of length steps+1
    (double[] path, double[]? antiPath) Generate(double S0, double r, double q, double sigma, double dt, int steps);
}

public sealed class GbmExactPathGen : IPathGen
{
    private readonly Random _rng;
    private readonly bool _antithetic;

    public GbmExactPathGen(int seed, bool antithetic = false)
    {
        _rng = new Random(seed);
        _antithetic = antithetic;
    }

    public (double[] path, double[]? antiPath) Generate(double S0, double r, double q, double sigma, double dt, int steps)
    {
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double vol = sigma * Math.Sqrt(dt);

        var p = new double[steps + 1]; p[0] = S0;
        double[]? a = _antithetic ? new double[steps + 1] : null;
        if (a != null) a[0] = S0;

        for (int i = 0; i < steps; i++)
        {
            double z = StdNormal();
            double zAnti = -z;

            p[i + 1] = p[i] * Math.Exp(drift + vol * z);
            if (a != null) a[i + 1] = a[i] * Math.Exp(drift + vol * zAnti);
        }
        return (p, a);

        double StdNormal()
        {
            // Box–Muller
            double u1 = 1.0 - _rng.NextDouble();
            double u2 = 1.0 - _rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
    }
}

public static class DeltaCv
{
    // payoffDelta(S, K, r, q, sigma, tau)
    public static double HedgingControlVariate(
        double[] path, double r, double q, double dt,
        Func<double, double, double, double, double, double, double> payoffDelta,
        double K, double sigma, double T)
    {
        double sum = 0.0;
        double dfDrift = Math.Exp((r - q) * dt);
        int m = path.Length - 1;

        for (int i = 0; i < m; i++)
        {
            double t = i * dt;
            double tau = Math.Max(T - t, 1e-12);
            double Si = path[i];
            double Si1 = path[i + 1];
            double delta = payoffDelta(Si, K, r, q, sigma, tau);

            // Zero-mean martingale increment under Q for S_t e^{-(r-q)t}
            sum += delta * (Si1 - dfDrift * Si);
        }
        return sum;
    }
}

// Simple struct for Welford accumulators we can merge safely.
internal struct RunningStats
{
    public long Count;
    public double Mean;
    public double M2;

    public void Push(double x)
    {
        Count++;
        double delta = x - Mean;
        Mean += delta / Count;
        M2 += delta * (x - Mean);
    }

    // Chan et al. merge of two Welford accumulators
    public static RunningStats Merge(RunningStats a, RunningStats b)
    {
        if (a.Count == 0) return b;
        if (b.Count == 0) return a;

        double delta = b.Mean - a.Mean;
        long n = a.Count + b.Count;

        RunningStats res = new RunningStats
        {
            Count = n,
            Mean = a.Mean + delta * b.Count / n,
            M2 = a.M2 + b.M2 + delta * delta * a.Count * b.Count / n
        };
        return res;
    }
}

public sealed class McPricer
{
    private readonly bool _useParallel;
    private readonly int _maxDegree;

    /// <param name="useParallel">Enable/disable multithreaded execution.</param>
    /// <param name="maxDegreeOfParallelism">Optional cap on threads; defaults to Environment.ProcessorCount when null or <=0.</param>
    public McPricer(bool useParallel = false, int? maxDegreeOfParallelism = null)
    {
        _useParallel = useParallel;
        _maxDegree = (maxDegreeOfParallelism.HasValue && maxDegreeOfParallelism.Value > 0)
            ? maxDegreeOfParallelism.Value
            : Environment.ProcessorCount;
    }

    // payoff(ST)
    public double PriceEuropean(
        double S0, double K, double r, double q, double sigma,
        double T, int steps, int nSamples, VarianceReduction vr,
        Func<double, double> payoff,
        Func<double, double, double, double, double, double, double> payoffDelta,
        out double stdError)
    {
        bool useAnti = vr.HasFlag(VarianceReduction.Antithetic);
        bool useCv = vr.HasFlag(VarianceReduction.DeltaControlVariate);

        if (useAnti && (nSamples % 2 != 0))
            throw new ArgumentException("nSamples must be even when Antithetic is enabled.");

        int nUnits = useAnti ? nSamples / 2 : nSamples; // each unit = one path or one antithetic pair
        double dt = T / steps;
        double disc = Math.Exp(-r * T);

        // SERIAL path (unchanged baseline)
        if (!_useParallel || nUnits <= 1)
        {
            var gen = new GbmExactPathGen(seed: 12345, antithetic: useAnti);
            double mean = 0.0, m2 = 0.0;
            int k = 0;

            for (int n = 0; n < nUnits; n++)
            {
                var (p, a) = gen.Generate(S0, r, q, sigma, dt, steps);

                double y1 = payoff(p[^1]);
                if (useCv) y1 -= DeltaCv.HedgingControlVariate(p, r, q, dt, payoffDelta, K, sigma, T);

                double y = y1;
                if (useAnti && a != null)
                {
                    double y2 = payoff(a[^1]);
                    if (useCv) y2 -= DeltaCv.HedgingControlVariate(a, r, q, dt, payoffDelta, K, sigma, T);
                    y = 0.5 * (y1 + y2);
                }

                k++;
                double deltaM = y - mean;
                mean += deltaM / k;
                m2 += deltaM * (y - mean);
            }

            double var = k > 1 ? m2 / (k - 1) : 0.0;
            double price = disc * mean;
            stdError = disc * Math.Sqrt(var / k);
            return price;
        }

        // PARALLEL path: only the unit loop is parallelized. RNG & payoff impls are unchanged.
        // Deterministic seeding: make each unit's generator from (baseSeed + unitIndex).
        const int baseSeed = 12345;

        var range = Partitioner.Create(0, nUnits, Math.Max(1, nUnits / (_maxDegree * 4)));
        var po = new ParallelOptions { MaxDegreeOfParallelism = _maxDegree };

        RunningStats globalStats = default;

        object mergeLock = new object();

        Parallel.ForEach(range, po, () => default(RunningStats), (slice, state, localStats) =>
        {
            for (int n = slice.Item1; n < slice.Item2; n++)
            {
                // One generator per unit for determinism across thread counts
                var gen = new GbmExactPathGen(seed: baseSeed + n, antithetic: useAnti);
                var (p, a) = gen.Generate(S0, r, q, sigma, dt, steps);

                double y1 = payoff(p[^1]);
                if (useCv) y1 -= DeltaCv.HedgingControlVariate(p, r, q, dt, payoffDelta, K, sigma, T);

                double y = y1;
                if (useAnti && a != null)
                {
                    double y2 = payoff(a[^1]);
                    if (useCv) y2 -= DeltaCv.HedgingControlVariate(a, r, q, dt, payoffDelta, K, sigma, T);
                    y = 0.5 * (y1 + y2);
                }

                localStats.Push(y);
            }
            return localStats;
        },
        localStats =>
        {
            // merge local into global
            lock (mergeLock)
            {
                globalStats = RunningStats.Merge(globalStats, localStats);
            }
        });

        long kTot = globalStats.Count;
        double varTot = kTot > 1 ? globalStats.M2 / (kTot - 1) : 0.0;
        double pricePar = disc * globalStats.Mean;
        stdError = disc * Math.Sqrt(varTot / kTot);
        return pricePar;
    }
}

public static class Program
{
    static void Main()
    {
        // Instrument and model parameters
        double S0 = 100.0;
        double K = 100.0;
        double r = 0.05;
        double q = 0.00;     // set dividend yield if needed
        double sigma = 0.2;
        double T = 1.0;      // in years

        int steps = 64;
        int nSamples = 200_000;

        // Payoff and delta for a European Call
        Func<double, double> payoffCall = ST => Math.Max(ST - K, 0.0);
        Func<double, double, double, double, double, double, double> deltaCall =
            (S, K_, r_, q_, sig_, tau_) => BlackScholes.CallDelta(S, K_, r_, q_, sig_, tau_);

        // Ground truth (closed-form)
        double bs = BlackScholes.PriceCall(S0, K, r, q, sigma, T);

        Console.WriteLine("European Call via Monte Carlo with Variance Reduction");
        Console.WriteLine($"S0={S0}, K={K}, r={r}, q={q}, sigma={sigma}, T={T}");
        Console.WriteLine($"Black–Scholes (ref): {bs:F6}");
        Console.WriteLine();

        // Baseline: single-threaded (for apples-to-apples comparison)
        var pricerSerial = new McPricer(useParallel: false);
        Report(pricerSerial, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.None,        "None (1 thread)");
        Report(pricerSerial, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.Antithetic,  "Antithetic (1 thread)");
        Report(pricerSerial, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.DeltaControlVariate, "Delta-CV (1 thread)");
        Report(pricerSerial, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.Antithetic | VarianceReduction.DeltaControlVariate, "Anti + Delta-CV (1 thread)");

        Console.WriteLine();

        // Parallel: auto-scales to Environment.ProcessorCount cores
        int cores = Environment.ProcessorCount;
        var pricerParallel = new McPricer(useParallel: true); // defaults to cores
        Report(pricerParallel, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.None,        $"None (parallel x{cores})");
        Report(pricerParallel, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.Antithetic,  $"Antithetic (parallel x{cores})");
        Report(pricerParallel, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.DeltaControlVariate, $"Delta-CV (parallel x{cores})");
        Report(pricerParallel, S0, K, r, q, sigma, T, steps, nSamples, payoffCall, deltaCall, VarianceReduction.Antithetic | VarianceReduction.DeltaControlVariate, $"Anti + Delta-CV (parallel x{cores})");
    }

    private static void Report(
        McPricer pricer,
        double S0, double K, double r, double q, double sigma, double T, int steps, int nSamples,
        Func<double, double> payoff,
        Func<double, double, double, double, double, double, double> payoffDelta,
        VarianceReduction vr, string label)
    {
        var t0 = DateTime.UtcNow;
        double price = pricer.PriceEuropean(S0, K, r, q, sigma, T, steps, nSamples, vr, payoff, payoffDelta, out double se);
        var t1 = DateTime.UtcNow;
        double ms = (t1 - t0).TotalMilliseconds;

        Console.WriteLine($"{label,-28}  Price = {price,10:F6}   StdErr = {se,10:F6}   Time = {ms,8:F1} ms");
    }
}
