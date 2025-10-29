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

public enum BarrierType { UpAndOut, DownAndOut, UpAndIn, DownAndIn }
public enum LookbackType { FixedStrike, FloatingStrike }

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
            sum += delta * (Si1 - dfDrift * Si);
        }
        return sum;
    }
}

// ===================== PAYOFFS (OOP) =====================
public interface IPathPayoff { double Payoff(double[] path); }

public abstract class OptionBase : IPathPayoff
{
    public double S0 { get; }
    public double K { get; }
    public double r { get; }
    public double q { get; }
    public double sigma { get; }
    public double T { get; }
    public int Steps { get; }

    protected OptionBase(double S0, double K, double r, double q, double sigma, double T, int steps)
    { S0 = S0; this.K = K; this.r = r; this.q = q; this.sigma = sigma; this.T = T; this.Steps = steps; }

    public abstract double Payoff(double[] path);
}

public sealed class EuropeanCall : IPathPayoff
{
    private readonly double _K;
    public EuropeanCall(double K) { _K = K; }
    public double Payoff(double[] path) => Math.Max(path[^1] - _K, 0.0);
}

public sealed class EuropeanPut : IPathPayoff
{
    private readonly double _K;
    public EuropeanPut(double K) { _K = K; }
    public double Payoff(double[] path) => Math.Max(_K - path[^1], 0.0);
}

public sealed class AsianArithmetic : IPathPayoff
{
    private readonly double _K; private readonly bool _isCall; private readonly bool _useAveragingFromZero;
    public AsianArithmetic(double K, bool isCall = true, bool useAveragingFromTimeZero = true)
    { _K = K; _isCall = isCall; _useAveragingFromZero = useAveragingFromTimeZero; }

    public double Payoff(double[] path)
    {
        int start = _useAveragingFromZero ? 0 : 1; // exclude S0 if desired
        double sum = 0.0; int n = 0;
        for (int i = start; i < path.Length; i++) { sum += path[i]; n++; }
        double avg = sum / n;
        double sign = _isCall ? 1.0 : -1.0;
        return Math.Max(sign * (avg - _K), 0.0);
    }
}

public sealed class DigitalOption : IPathPayoff
{
    private readonly double _K; private readonly bool _isCall; private readonly double _cash;
    public DigitalOption(double K, bool isCall = true, double cashPayout = 1.0)
    { _K = K; _isCall = isCall; _cash = cashPayout; }
    public double Payoff(double[] path)
    {
        double ST = path[^1];
        bool hit = _isCall ? (ST > _K) : (ST < _K);
        return hit ? _cash : 0.0;
    }
}

public sealed class BarrierOption : IPathPayoff
{
    private readonly double _K; private readonly bool _isCall; private readonly BarrierType _type; private readonly double _barrier; private readonly double _rebate;
    public BarrierOption(double K, bool isCall, BarrierType type, double barrier, double rebate = 0.0)
    { _K = K; _isCall = isCall; _type = type; _barrier = barrier; _rebate = rebate; }

    public double Payoff(double[] path)
    {
        bool touched = false;
        for (int i = 0; i < path.Length; i++)
        {
            double S = path[i];
            if ((_type == BarrierType.UpAndOut || _type == BarrierType.UpAndIn) && S >= _barrier) { touched = true; break; }
            if ((_type == BarrierType.DownAndOut || _type == BarrierType.DownAndIn) && S <= _barrier) { touched = true; break; }
        }

        double intrinsic = _isCall ? Math.Max(path[^1] - _K, 0.0) : Math.Max(_K - path[^1], 0.0);

        switch (_type)
        {
            case BarrierType.UpAndOut:
            case BarrierType.DownAndOut:
                return touched ? _rebate : intrinsic;
            case BarrierType.UpAndIn:
            case BarrierType.DownAndIn:
                return touched ? intrinsic : _rebate;
            default: return 0.0;
        }
    }
}

public sealed class LookbackOption : IPathPayoff
{
    private readonly bool _isCall; private readonly LookbackType _type; private readonly double _K;
    public LookbackOption(bool isCall, LookbackType type, double K = 0.0) { _isCall = isCall; _type = type; _K = K; }

    public double Payoff(double[] path)
    {
        double Smin = double.MaxValue, Smax = double.MinValue;
        foreach (double S in path) { if (S < Smin) Smin = S; if (S > Smax) Smax = S; }
        double ST = path[^1];

        if (_type == LookbackType.FloatingStrike)
        {
            return _isCall ? Math.Max(ST - Smin, 0.0) : Math.Max(Smax - ST, 0.0);
        }
        else // Fixed strike
        {
            return _isCall ? Math.Max(Smax - _K, 0.0) : Math.Max(_K - Smin, 0.0);
        }
    }
}

public sealed class RangeOption : IPathPayoff
{
    // Payoff = max( (S_max - S_min) - K, 0 ), set K=0 for pure range
    private readonly double _K;
    public RangeOption(double K = 0.0) { _K = K; }
    public double Payoff(double[] path)
    {
        double Smin = double.MaxValue, Smax = double.MinValue;
        foreach (double S in path) { if (S < Smin) Smin = S; if (S > Smax) Smax = S; }
        return Math.Max((Smax - Smin) - _K, 0.0);
    }
}

// ===================== STATS =====================
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

// ===================== MONTE CARLO PRICER =====================
public sealed class McPricer
{
    private readonly bool _useParallel;
    private readonly int _maxDegree;
    private readonly int _baseSeed;

    public McPricer(bool useParallel = false, int? maxDegreeOfParallelism = null, int baseSeed = 12345)
    {
        _useParallel = useParallel;
        _maxDegree = (maxDegreeOfParallelism.HasValue && maxDegreeOfParallelism.Value > 0)
            ? maxDegreeOfParallelism.Value
            : Environment.ProcessorCount;
        _baseSeed = baseSeed;
    }

    public double Price(
        double S0, double r, double q, double sigma, double T,
        int steps, int nSamples, VarianceReduction vr,
        IPathPayoff payoff,
        Func<double, double, double, double, double, double, double>? payoffDelta,
        double KforCV,
        out double stdError)
    {
        bool useAnti = vr.HasFlag(VarianceReduction.Antithetic);
        bool useCv = vr.HasFlag(VarianceReduction.DeltaControlVariate);
        if (useAnti && (nSamples % 2 != 0)) throw new ArgumentException("nSamples must be even when Antithetic is enabled.");

        int nUnits = useAnti ? nSamples / 2 : nSamples; // each unit = one path or one antithetic pair
        double dt = T / steps;
        double disc = Math.Exp(-r * T);

        if (!_useParallel || nUnits <= 1)
        {
            var gen = new GbmExactPathGen(seed: _baseSeed, antithetic: useAnti);
            double mean = 0.0, m2 = 0.0; int k = 0;
            for (int n = 0; n < nUnits; n++)
            {
                var (p, a) = gen.Generate(S0, r, q, sigma, dt, steps);
                double y1 = payoff.Payoff(p);
                if (useCv && payoffDelta != null)
                    y1 -= DeltaCv.HedgingControlVariate(p, r, q, dt, payoffDelta, KforCV, sigma, T);

                double y = y1;
                if (useAnti && a != null)
                {
                    double y2 = payoff.Payoff(a);
                    if (useCv && payoffDelta != null)
                        y2 -= DeltaCv.HedgingControlVariate(a, r, q, dt, payoffDelta, KforCV, sigma, T);
                    y = 0.5 * (y1 + y2);
                }

                k++;
                double deltaM = y - mean;
                mean += deltaM / k;
                m2 += deltaM * (y - mean);
            }
            double var = k > 1 ? m2 / (k - 1) : 0.0;
            stdError = disc * Math.Sqrt(var / k);
            return disc * mean;
        }
        else
        {
            var range = Partitioner.Create(0, nUnits, Math.Max(1, nUnits / (_maxDegree * 4)));
            var po = new ParallelOptions { MaxDegreeOfParallelism = _maxDegree };
            RunningStats globalStats = default; object mergeLock = new object();

            Parallel.ForEach(range, po, () => default(RunningStats), (slice, state, localStats) =>
            {
                for (int n = slice.Item1; n < slice.Item2; n++)
                {
                    var gen = new GbmExactPathGen(seed: _baseSeed + n, antithetic: useAnti);
                    var (p, a) = gen.Generate(S0, r, q, sigma, dt, steps);
                    double y1 = payoff.Payoff(p);
                    if (useCv && payoffDelta != null)
                        y1 -= DeltaCv.HedgingControlVariate(p, r, q, dt, payoffDelta, KforCV, sigma, T);

                    double y = y1;
                    if (useAnti && a != null)
                    {
                        double y2 = payoff.Payoff(a);
                        if (useCv && payoffDelta != null)
                            y2 -= DeltaCv.HedgingControlVariate(a, r, q, dt, payoffDelta, KforCV, sigma, T);
                        y = 0.5 * (y1 + y2);
                    }
                    localStats.Push(y);
                }
                return localStats;
            },
            localStats => { lock (mergeLock) { globalStats = RunningStats.Merge(globalStats, localStats); } });

            long kTot = globalStats.Count;
            double varTot = kTot > 1 ? globalStats.M2 / (kTot - 1) : 0.0;
            stdError = Math.Exp(-r * T) * Math.Sqrt(varTot / kTot);
            return Math.Exp(-r * T) * globalStats.Mean;
        }
    }

    // ======= Greek estimators via CRN bump-and-revalue (central differences) =======
    public (double Delta, double Gamma, double Vega, double Rho, double Theta) Greeks(
        double S0, double r, double q, double sigma, double T,
        int steps, int nSamples, VarianceReduction vr,
        IPathPayoff payoff,
        Func<double, double, double, double, double, double, double>? payoffDeltaForCV,
        double KforCV,
        double dS = 1e-2, double dSigma = 1e-4, double dR = 1e-4, double dT = 1.0 / 365.0)
    {
        // Use same baseSeed to enforce CRN
        double Price0 = Price(S0, r, q, sigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);
        double PriceSUp = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0 + dS, r, q, sigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);
        double PriceSDown = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0 - dS, r, q, sigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);

        double PriceVUp = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0, r, q, sigma + dSigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);
        double PriceVDown = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0, r, q, sigma - dSigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);

        double PriceRUp = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0, r + dR, q, sigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);
        double PriceRDown = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0, r - dR, q, sigma, T, steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);

        double PriceTUp = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0, r, q, sigma, Math.Max(T + dT, 1e-9), steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);
        double PriceTDown = new McPricer(_useParallel, _maxDegree, _baseSeed).Price(S0, r, q, sigma, Math.Max(T - dT, 1e-9), steps, nSamples, vr, payoff, payoffDeltaForCV, KforCV, out _);

        double delta = (PriceSUp - PriceSDown) / (2.0 * dS);
        double gamma = (PriceSUp - 2.0 * Price0 + PriceSDown) / (dS * dS);
        double vega = (PriceVUp - PriceVDown) / (2.0 * dSigma);
        double rho  = (PriceRUp - PriceRDown) / (2.0 * dR);
        double theta = (PriceTDown - PriceTUp) / (2.0 * dT); // dPrice/dT, negative is usual
        return (delta, gamma, vega, rho, theta);
    }
}

public static class Demo
{
    public static void Main()
    {
        double S0 = 100.0, K = 100.0, r = 0.03, q = 0.00, sigma = 0.2, T = 1.0;
        int steps = 252, nSamples = 200_000;
        var pricer = new McPricer(useParallel: true);

        // Control variate uses European call delta by default
        Func<double, double, double, double, double, double, double> cvDelta = (S, K_, r_, q_, sig_, tau_) => BlackScholes.CallDelta(S, K_, r_, q_, sig_, tau_);

        // 1) Asian arithmetic average CALL
        var asian = new AsianArithmetic(K: K, isCall: true);
        double pAsian = pricer.Price(S0, r, q, sigma, T, steps, nSamples,
            VarianceReduction.Antithetic | VarianceReduction.DeltaControlVariate,
            asian, cvDelta, K, out double seAsian);

        var gAsian = pricer.Greeks(S0, r, q, sigma, T, steps, nSamples,
            VarianceReduction.Antithetic | VarianceReduction.DeltaControlVariate,
            asian, cvDelta, K);

        Console.WriteLine($"Asian (arith) Call: Price={pAsian:F6}, SE={seAsian:F6}, Delta={gAsian.Delta:F5}, Gamma={gAsian.Gamma:F5}, Vega={gAsian.Vega:F5}, Rho={gAsian.Rho:F5}, Theta={gAsian.Theta:F5}");

        // 2) Digital cash-or-nothing CALL (cash=1)
        var digital = new DigitalOption(K: K, isCall: true, cashPayout: 1.0);
        double pDigital = pricer.Price(S0, r, q, sigma, T, steps, nSamples,
            VarianceReduction.Antithetic, digital, cvDelta, K, out double seDigital);
        var gDigital = pricer.Greeks(S0, r, q, sigma, T, steps, nSamples, VarianceReduction.Antithetic, digital, cvDelta, K);
        Console.WriteLine($"Digital Call: Price={pDigital:F6}, SE={seDigital:F6}, Delta={gDigital.Delta:F5}, Gamma={gDigital.Gamma:F5}, Vega={gDigital.Vega:F5}, Rho={gDigital.Rho:F5}, Theta={gDigital.Theta:F5}");

        // 3) Barrier: Up-and-Out Call with barrier=120, no rebate
        var barrier = new BarrierOption(K: K, isCall: true, type: BarrierType.UpAndOut, barrier: 120.0);
        double pBarrier = pricer.Price(S0, r, q, sigma, T, steps, nSamples,
            VarianceReduction.Antithetic, barrier, cvDelta, K, out double seBarrier);
        var gBarrier = pricer.Greeks(S0, r, q, sigma, T, steps, nSamples, VarianceReduction.Antithetic, barrier, cvDelta, K);
        Console.WriteLine($"Barrier UO Call (B=120): Price={pBarrier:F6}, SE={seBarrier:F6}, Delta={gBarrier.Delta:F5}, Gamma={gBarrier.Gamma:F5}, Vega={gBarrier.Vega:F5}, Rho={gBarrier.Rho:F5}, Theta={gBarrier.Theta:F5}");

        // 4) Lookback Floating-Strike Call
        var lookback = new LookbackOption(isCall: true, type: LookbackType.FloatingStrike);
        double pLb = pricer.Price(S0, r, q, sigma, T, steps, nSamples,
            VarianceReduction.Antithetic, lookback, cvDelta, K, out double seLb);
        var gLb = pricer.Greeks(S0, r, q, sigma, T, steps, nSamples, VarianceReduction.Antithetic, lookback, cvDelta, K);
        Console.WriteLine($"Lookback Floating Call: Price={pLb:F6}, SE={seLb:F6}, Delta={gLb.Delta:F5}, Gamma={gLb.Gamma:F5}, Vega={gLb.Vega:F5}, Rho={gLb.Rho:F5}, Theta={gLb.Theta:F5}");

        // 5) Range option (K=0 => pure range)
        var range = new RangeOption(K: 0.0);
        double pRange = pricer.Price(S0, r, q, sigma, T, steps, nSamples,
            VarianceReduction.Antithetic, range, cvDelta, K, out double seRange);
        var gRange = pricer.Greeks(S0, r, q, sigma, T, steps, nSamples, VarianceReduction.Antithetic, range, cvDelta, K);
        Console.WriteLine($"Range (pure): Price={pRange:F6}, SE={seRange:F6}, Delta={gRange.Delta:F5}, Gamma={gRange.Gamma:F5}, Vega={gRange.Vega:F5}, Rho={gRange.Rho:F5}, Theta={gRange.Theta:F5}");
    }
}
