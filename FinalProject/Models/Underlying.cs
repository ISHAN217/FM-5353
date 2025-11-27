using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;

namespace FinalProject.Models
{
    public class Underlying
    {
        [Key]
        public int Id { get; set; }

        [Required]
        public string Symbol { get; set; } = string.Empty; // e.g., AAPL

        [Required]
        public string CompanyName { get; set; } = string.Empty;

        // Current spot price of the stock
        // We use decimal for money/financial data to avoid rounding errors
        [Column(TypeName = "decimal(18,4)")]
        public decimal CurrentPrice { get; set; }

        // Volatility (Sigma) - crucial for Monte Carlo
        [Column(TypeName = "decimal(10,4)")]
        public decimal Volatility { get; set; }

        // This allows us to see all derivatives attached to this stock
        public List<Derivative> Derivatives { get; set; } = new();
    }
}