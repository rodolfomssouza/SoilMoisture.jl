using SoilMoisture
using Plots

# %% Parameters --------------------------------------------------------------
s1 = 0.10
s2 = 0.20
s3 = 0.45
s4 = 1.00
sh = 0.08
sw = 0.11
sstar = 0.33
sfc = 0.35
emax = 0.50
ew = emax / 10
b = (4.05 * 2) + 3
ks = 201.0
n = 0.37
zr = 40.0
dt = 1 / 48
rain1 = 0.0
rain2 = 10.0
rain3 = zeros(10)

params = (
    sh = sh,
    sw = sw,
    sstar = sstar,
    sfc = sfc,
    emax = emax,
    ew = ew,
    b = b,
    ks = ks,
    n = n,
    zr = zr,
    dt = dt,
)

# %% Rainfall series ---------------------------------------------------------
years = 101
days = years * 365
number_events = Int(days * 1 / dt)
α = 0.90
λ = 0.35
rain = rainfall_poisson(number_events, α, λ * dt)


# Run model ------------------------------------------------------------------
s5 = 0.25
res1 = solve_swb(rain, s5, params)
res2 = dt2daily(res1)

# %% Calculate PR and ET -----------------------------------------------------
bd = 1.68
a = -5.76
b = 5.63
c = -15.32
rp = @. soil_rp(res2.s_mean * n, bd, a, b, c)


# %% Plot histogram of PR ----------------------------------------------------
p1 = histogram(
    res2.s_mean,
    bins = 100,
    label = "Soil moisture",
    xlabel = "Soil moisture",
    ylabel = "Density",
    normalized = :true,
)

p2 = histogram(
    rp,
    bins = 100,
    label = "PR",
    xlabel = "PR (MPa)",
    ylabel = "Density",
    normalized = :true,
)

plot(p1, p2, layout = (2, 1), size = (800, 600))

