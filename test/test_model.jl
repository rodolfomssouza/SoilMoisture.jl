using SoilMoisture

# %% Parameters --------------------------------------------------------------
s1 = 0.10
s2 = 0.20
s3 = 0.45
s4 = 1.00
sh = 0.12
sw = 0.25
sstar = 0.50
sfc = 0.65
emax = 0.50
ew = emax / 10
b = 4.0
ks = 20
n = 0.50
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
    dt = dt
         )

rain = zeros(48 * 10)
s5 = 0.70

res1 = solve_swb(rain, s5, params)
res2 = dt2daily(res1)

