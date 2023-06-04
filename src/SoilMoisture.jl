module SoilMoisture

# Write your package code here.
# %% Packages ----------------------------------------------------------------
using DataFrames
using Distributions


# %% Exports -----------------------------------------------------------------
export rainfall_poisson
export evapotranspiration
export leakage
export water_loss
export soil_water_balance
export solve_swb
export dt2daily
export soil_rp 


# %% Functions ---------------------------------------------------------------

"""
`rainfall_poisson(n, α, λ)`

Generates rainfall using Poisson process.

# Arguments
- `n::Integer`: number of events.
- `α::Float64`: mean rain depth.
- `λ::Float64`: rain event interval (1/d).

If rainfall needs to be similated in timescale smaller than daily, multiply `λ` by dt.

# Example
```
nevents = 30
α = 0.75
λ = 0.45
rain = rainfall_poisson(n, α, λ)
```
"""
function rainfall_poisson(nevents::Int64, α::Float64, λ::Float64)
    rain::Vector{Float64} = zeros(nevents)
    r1::Float64 = 0.0
    r2::Float64 = 0.0
    for i = 1:nevents
        r1 = rand(Uniform(0, 1))
        if r1 < λ
            r2 = rand(Uniform(0, 1))
            rain[i] = -α * log((1 - r2))
        else
            rain[i] = 0
        end
    end
    return rain
end



"""
`evapotranspirationi(s, sh, sw, sstar, emax, ew)`

Computes evapotranspiration losses

# Arguments
- `s`: soil moisture
- `sh`: soil moisture at hygroscopic point.
- `sw`: soil moisture at wilting point.
- `sstar`: soil moisture below field capacity.
- `emax`: maximum evapotranspiration rate.
- `ew`: soil evaporation rate

# Example
```
s = 0.45
sh = 0.12
sw = 0.20
sstar = 0.55
emax = 0.50
ew = 0.05

et = evapotranspiration(s, sh, sw, sstar, emax, ew)
```
"""
function evapotranspiration(s, sh, sw, sstar, emax, ew)
    if s <= sh
        return 0.0
    elseif s <= sw
        return ew * (s - sh) / (sw - sh)
    elseif s <= sstar
        return ew + (emax - ew) * (s - sw) / (sstar - sw)
    else
        return emax
    end
end


"""
`leakage(s, sfc, b, ks)`

Computes the soil leakage

# Arguments
- `s`: soil moisture.
- `sfc`: soil moisture at field capacity.
- `b`: leakage curve exponent.
- `ks`: hydralic conductivity.

# Example
```
s = 0.60
sfc = 0.55
b = 4.5
ks = 150

lk = leakage(s, sfc, b, ks)
```
"""
function leakage(s, sfc, b, ks)::Float64
    if s <= sfc
        lk = 0
    else
        lk = ks * s^b
    end
    return lk
end


"""
`water_loss(s, sh, sw, sstar, sfc, emax, ew, b, ks)`

Water loss function

Combines evapotranspiration and leakage as a function of soil moisture.

# Arguments
- `s`: soil moisture.
- `sh`: soil moisture at hygroscopic point.
- `sw`: soil moisture at wilting point.
- `sstar`: soil moisture below field capacity.
- `sfc`: soil moisture at field capacity.
- `emax`: maximum evapotranspiration rate.
- `ew`: soil evaporation rate.
- `b`: leakage curve exponent.
- `ks`: hydralic conductivity.
"""
function water_loss(s, sh, sw, sstar, sfc, emax, ew, b, ks)::Float64
    et = evapotranspiration(s, sh, sw, sstar, emax, ew)
    lk = leakage(s, sfc, b, ks)
    return et + lk
end


"""
Computes the soil water balance

This function assumes that the canopy interception is negligible.

# Arguments
- `rain`: rainfall.
- `s`: soil moisture.
- `sh`: soil moisture at hygroscopic point.
- `sw`: soil moisture at wilting point.
- `sstar`: soil moisture below field capacity.
- `sfc`: soil moisture at field capacity.
- `b`: leakage curve exponent.
- `ks`: hydralic conductivity.
- `n`: porosity.
- `zr`: root zone depth.
- `emax`: maximum evapotranspiration rate.
- `ew`: soil evaporation rate.
- `dt`: time step.
"""
function soil_water_balance(
    rain,
    s,
    sh,
    sw,
    sstar,
    sfc,
    b,
    ks,
    n,
    zr,
    emax,
    ew,
    dt,
)
    # Soil water storage
    nzr = n * zr

    # Add rainfall
    s += rain / nzr

    # Check for runoff
    runoff = max(s - 1.0, 0.0)
    s = min(s, 1.0)

    # Evapotranspiration
    et = evapotranspiration(s, sh, sw, sstar, emax, ew) * dt
    s -= et / nzr

    # Leakage
    lk = leakage(s, sfc, b, ks) * dt
    s -= lk / nzr

    # Return fluxes
    res = [s, et, lk, runoff]
    return res
end


"""
Solve soil water balance model

# Arguments
- `rain`: a vector of rainfall series at the same time step the model will be solved.
- `s`: initial guess for the soil moisture value.
- `params`: a dictionary with all parameters needed to for `soil_water_balance()` function.
"""
function solve_swb(rain, s, params)

    # Unzip parameters
    sh = params[:sh]
    sw = params[:sw]
    sstar = params[:sstar]
    sfc = params[:sfc]
    b = params[:b]
    ks = params[:ks]
    n = params[:n]
    zr = params[:zr]
    emax = params[:emax]
    ew = params[:ew]
    dt = params[:dt]

    # Create object to receive results
    nr = length(rain) + 1
    s_res = zeros(nr)
    et_res = zeros(nr)
    lk_res = zeros(nr)
    runoff_res = zeros(nr)
    s_res[1] = s

    # Create vector for time
    ndays = Int((nr - 1) * dt)
    days = [1:ndays;]
    days = repeat(days, inner = Int(1 / dt))
    days = vcat(0, days)
    days_cont = @. [(1 + dt):dt:(ndays + 1);] - dt

    # Solve model
    for i in eachindex(rain)
        sol = soil_water_balance(
            rain[i],
            s_res[i],
            sh,
            sw,
            sstar,
            sfc,
            b,
            ks,
            n,
            zr,
            emax,
            ew,
            dt,
        )
        s_res[i + 1] = sol[1]
        et_res[i + 1] = sol[2]
        lk_res[i + 1] = sol[3]
        runoff_res[i + 1] = sol[4]
    end

    # Create dataframe to export results
    df = DataFrame(
        Days = days,
        DaysCont = vcat(0, days_cont),
        Rain = vcat(0, rain),
        Q = runoff_res,
        Lk = lk_res,
        ET = et_res,
        s = s_res,
    )
    return df
end


"""
Convert soil water balance results to daily time-scale
"""
function dt2daily(df)
    df1 = combine(
        groupby(df, :Days),
        :Rain => sum,
        :Q => sum,
        :s => mean,
        :Lk => sum,
        :ET => sum,
        skipmissing = true,
    )
    return df1
end


"""
Soil penetration resistance

`soil_rp(s, bulk_density, a, b, c)`

# Arguments
- `s`: soil moisture.
- `bulk_density`: soil bulk density.
- `a`, `b`, and `c`: soil penetration resistance parameters.

# Example
```
s = 0.45
bd = 1.68
a = -5.75
b = 6.45
c = -15.30
rp = soil_rp(s, bulk_density, a, b, c)
```
"""
function soil_rp(s, bulk_density, a, b, c)
    rp = exp(a + b * bulk_density + c * s)
    return rp
end

end
