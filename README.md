🚀 Full-Stack Monte Carlo Option Pricer

A high-performance financial application that calculates Market Value, Delta, Gamma, Vega, Theta, and Rho for European Options using Monte Carlo simulations.

🏗 Architecture

This project uses a modern Client-Server architecture:

Backend: ASP.NET Core 9 Web API (C#)

Database: PostgreSQL (Entity Framework Core)

Math Engine: Custom Monte Carlo Simulation Service (10,000 iterations)

Frontend: Blazor WebAssembly (Real-time Dashboard)

✨ Features

Real-Time Valuation: Calculates option prices on-the-fly using stochastic calculus.

Risk Analysis: Visualizes Delta (Directional Risk) and Gamma (Curvature) with dynamic CSS charts.

Advanced Greeks: Includes calculations for Vega, Theta, and Rho.

Trade Management: Supports both Call and Put options.

Persistent Storage: All trades are saved to a relational database.

🚀 How to Run

Prerequisites

.NET 9.0 SDK: Download Here

PostgreSQL: Download Here (or via Homebrew/Postgres.app)

Step 1: Start the Database & API (The Backend)

Open your terminal.

Navigate to the backend folder:

cd FinalProject



Run the application:

dotnet run



Look for the success message: Now listening on: http://localhost:5225

Note: Keep this terminal open! The backend must be running for the math to work.

Step 2: Start the Dashboard (The Frontend)

Open a new terminal window.

Navigate to the frontend folder:

cd Frontend



Run the website:

dotnet run



Look for the success message: Now listening on: http://localhost:5246

Step 3: Use the App

Open your web browser (Chrome, Safari, Edge).

Go to the Frontend URL: http://localhost:5246

Create a Trade:

Enter a ticker (e.g., "NVDA") and volatility (e.g., "0.45").

Select Call (Bet Up) or Put (Bet Down).

Click 🚀 Execute Trade.

Analyze Results:

The table will automatically update with the calculated Price and Greeks.

Look for the Risk Bar in the Delta column (Blue for Calls, Red for Puts).

🧮 The Math Behind It

The engine simulates 10,000 price paths using Geometric Brownian Motion:

$$ dS_t = \mu S_t dt + \sigma S_t dW_t $$

It calculates Greeks using the Finite Difference Method by "bumping" variables and re-running the simulation:

Delta ($\Delta$): Sensitivity to underlying price ($S \pm 1\%$)

Vega ($\nu$): Sensitivity to volatility ($\sigma + 1\%$)

Theta ($\Theta$): Sensitivity to time decay ($T - 1$ day)

Rho ($\rho$): Sensitivity to interest rates ($r + 1\%$)

🛠 Troubleshooting

"Connection Refused": Ensure PostgreSQL is running (brew services start postgresql).

"Failed to fetch": Ensure the API terminal is still running on port 5225.
