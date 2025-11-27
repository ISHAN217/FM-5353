using Microsoft.EntityFrameworkCore;
using FinalProject.Models;

namespace FinalProject.Data
{
    public class FinancialDbContext : DbContext
    {
        public FinancialDbContext(DbContextOptions<FinancialDbContext> options) : base(options) { }

        // These three lines tell Entity Framework to create these three tables
        public DbSet<Underlying> Underlyings { get; set; }
        public DbSet<Derivative> Derivatives { get; set; }
        public DbSet<Trade> Trades { get; set; }

        protected override void OnModelCreating(ModelBuilder modelBuilder)
        {
            // Optional: Seed some initial data so the app isn't empty when you start
            modelBuilder.Entity<Underlying>().HasData(
                new Underlying { 
                    Id = 1, 
                    Symbol = "AAPL", 
                    CompanyName = "Apple Inc.", 
                    CurrentPrice = 150.00m, 
                    Volatility = 0.25m 
                },
                new Underlying { 
                    Id = 2, 
                    Symbol = "TSLA", 
                    CompanyName = "Tesla Inc.", 
                    CurrentPrice = 200.00m, 
                    Volatility = 0.45m 
                }
            );
        }
    }
}