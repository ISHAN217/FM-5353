using Microsoft.AspNetCore.Mvc;
using Microsoft.EntityFrameworkCore;
using FinalProject.Data;
using FinalProject.Models;
using FinalProject.Services;

namespace FinalProject.Controllers
{
    [ApiController]
    [Route("api/[controller]")]
    public class PortfolioController : ControllerBase
    {
        private readonly FinancialDbContext _context;
        private readonly IPricingService _pricingService;

        public PortfolioController(FinancialDbContext context, IPricingService pricingService)
        {
            _context = context;
            _pricingService = pricingService;
        }

        [HttpGet("trades")]
        public async Task<ActionResult<IEnumerable<Trade>>> GetTrades()
        {
            return await _context.Trades
                .Include(t => t.Derivative)
                .ThenInclude(d => d.Underlying)
                .ToListAsync();
        }

        [HttpPost("trade")]
        public async Task<ActionResult<Trade>> CreateTrade(Trade trade)
        {
            _context.Trades.Add(trade);
            await _context.SaveChangesAsync();
            return CreatedAtAction(nameof(GetTrades), new { id = trade.Id }, trade);
        }

        [HttpPost("evaluate")]
        public async Task<ActionResult<List<ValuationResult>>> EvaluatePortfolio()
        {
            var trades = await _context.Trades
                .Include(t => t.Derivative)
                .ThenInclude(d => d.Underlying)
                .ToListAsync();

            var results = new List<ValuationResult>();

            foreach (var trade in trades)
            {
                var sim = _pricingService.Price(trade);

                results.Add(new ValuationResult
                {
                    TradeId = trade.Id,
                    Symbol = trade.Derivative?.Underlying?.Symbol ?? "Unknown",
                    Quantity = trade.Quantity,
                    MarketValue = sim.Price * trade.Quantity,
                    // Map all the Greeks
                    Delta = sim.Delta,
                    Gamma = sim.Gamma,
                    Vega = sim.Vega,
                    Theta = sim.Theta,
                    Rho = sim.Rho
                });
            }

            return Ok(results);
        }
    }

    public class ValuationResult
    {
        public int TradeId { get; set; }
        public string Symbol { get; set; } = string.Empty;
        public int Quantity { get; set; }
        public decimal MarketValue { get; set; }
        public decimal Delta { get; set; }
        public decimal Gamma { get; set; }
        public decimal Vega { get; set; }
        public decimal Theta { get; set; }
        public decimal Rho { get; set; }
    }
}