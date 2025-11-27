using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;

namespace FinalProject.Models
{
    public class Trade
    {
        [Key]
        public int Id { get; set; }

        [Required]
        public int Quantity { get; set; } // Positive = Buy, Negative = Sell

        public DateTime TradeDate { get; set; } = DateTime.UtcNow;

        [Required]
        public string Status { get; set; } = "Active";

        // Foreign Key: Links this trade to a specific Derivative
        public int DerivativeId { get; set; }
        public Derivative? Derivative { get; set; }
    }
}