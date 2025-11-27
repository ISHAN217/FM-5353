using FinalProject.Models;

namespace FinalProject.Services
{
    public class MonteCarloService : IPricingService
    {
        private const int SIMULATIONS = 10000;
        private readonly Random _rand = new Random();

        public SimulationResult Price(Trade trade)
        {
            var d = trade.Derivative;
            var u = d.Underlying;

            // 1. Base Variables
            double S = (double)u.CurrentPrice;
            double K = (double)d.StrikePrice;
            double r = (double)d.RiskFreeRate;
            double sigma = (double)u.Volatility;
            double T = (d.ExpirationDate - DateTime.Now).TotalDays / 365.0;
            if (T < 0) T = 0;
            bool isPut = d.Type.Contains("Put", StringComparison.OrdinalIgnoreCase);

            // 2. Calculate Base Price
            decimal price = RunSimulation(S, K, r, sigma, T, isPut);

            // 3. Calculate Greeks using Finite Difference Method (Bumping)
            
            // DELTA & GAMMA (Bump Price by 1%)
            double bumpS = S * 0.01;
            decimal priceUp = RunSimulation(S + bumpS, K, r, sigma, T, isPut);
            decimal priceDown = RunSimulation(S - bumpS, K, r, sigma, T, isPut);
            decimal delta = (priceUp - priceDown) / (2 * (decimal)bumpS);
            decimal gamma = (priceUp - 2 * price + priceDown) / (decimal)Math.Pow(bumpS, 2);

            // VEGA (Bump Volatility by 1%)
            double bumpSigma = 0.01;
            decimal priceVega = RunSimulation(S, K, r, sigma + bumpSigma, T, isPut);
            decimal vega = (priceVega - price) / (decimal)(bumpSigma * 100); // Scaled

            // THETA (Decrease Time by 1 day)
            double bumpT = 1.0 / 365.0;
            decimal priceTheta = RunSimulation(S, K, r, sigma, T - bumpT, isPut);
            decimal theta = (priceTheta - price); // Daily decay

            // RHO (Bump Rate by 1%)
            double bumpR = 0.01;
            decimal priceRho = RunSimulation(S, K, r + bumpR, sigma, T, isPut);
            decimal rho = (priceRho - price) / (decimal)(bumpR * 100);

            return new SimulationResult
            {
                Price = price,
                Delta = delta,
                Gamma = gamma,
                Vega = vega,
                Theta = theta,
                Rho = rho
            };
        }

        // Helper Method to run the loop
        private decimal RunSimulation(double S, double K, double r, double sigma, double T, bool isPut)
        {
            if (T <= 0) return isPut ? (decimal)Math.Max(K - S, 0) : (decimal)Math.Max(S - K, 0);

            double sumPayoff = 0;
            
            // Pre-calculate constants for speed
            double drift = (r - 0.5 * Math.Pow(sigma, 2)) * T;
            double volT = sigma * Math.Sqrt(T);

            for (int i = 0; i < SIMULATIONS; i++)
            {
                // Box-Muller Normal Distribution
                double u1 = 1.0 - _rand.NextDouble();
                double u2 = 1.0 - _rand.NextDouble();
                double z = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);

                double ST = S * Math.Exp(drift + volT * z);
                double payoff = isPut ? Math.Max(K - ST, 0) : Math.Max(ST - K, 0);
                sumPayoff += payoff;
            }

            double discountFactor = Math.Exp(-r * T);
            return (decimal)((sumPayoff / SIMULATIONS) * discountFactor);
        }
    }
}