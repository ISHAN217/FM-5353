using System.Text.Json.Serialization;
using Microsoft.AspNetCore.Builder;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;

var builder = WebApplication.CreateBuilder(args);

// Controllers + JSON enum as strings ("Call", "Put", "AsianArithmetic", etc.)
builder.Services
    .AddControllers()
    .AddJsonOptions(o =>
    {
        o.JsonSerializerOptions.Converters.Add(new JsonStringEnumConverter());
    });

// Relaxed CORS for localhost/Postman
builder.Services.AddCors(x => x.AddPolicy("permissive",
    b => b.AllowAnyOrigin()
          .AllowAnyMethod()
          .AllowAnyHeader()));

var app = builder.Build();

// HTTPS redirection disabled for local dev (per lecture)
 // app.UseHttpsRedirection();

app.UseCors("permissive");

app.MapControllers();

app.Run();

