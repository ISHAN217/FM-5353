using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;

namespace FinalProject.Models
{
    public class Derivative
    {
        [Key]
        public int Id { get; set; }

        [Required]
        public string Type { get; set; } = "European Call"; // Call, Put, etc.

        [Column(TypeName = "decimal(18,4)")]
        public decimal StrikePrice { get; set; }

        [Required]
        public DateTime ExpirationDate { get; set; }

        // Risk Free Rate (r)
        [Column(TypeName = "decimal(10,4)")]
        public decimal RiskFreeRate { get; set; }

        // Foreign Key: Links this derivative to a specific Underlying (stock)
        public int UnderlyingId { get; set; }
        public Underlying? Underlying { get; set; }

        public List<Trade> Trades { get; set; } = new();
    }
}