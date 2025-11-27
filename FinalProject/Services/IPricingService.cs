using FinalProject.Models;

namespace FinalProject.Services
{
    public interface IPricingService
    {
        SimulationResult Price(Trade trade);
    }

    public class SimulationResult
    {
        public decimal Price { get; set; }
        public decimal Delta { get; set; } // Sensitivity to Price
        public decimal Gamma { get; set; } // Curvature
        public decimal Vega { get; set; }  // Sensitivity to Volatility
        public decimal Theta { get; set; } // Time Decay
        public decimal Rho { get; set; }   // Sensitivity to Interest Rates
    }
}