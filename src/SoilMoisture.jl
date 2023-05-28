module SoilMoisture

# Write your package code here.
# %% Packages ----------------------------------------------------------------
using Distributions


# %% Exports -----------------------------------------------------------------
export rainfall_poisson
export evapotranspiration
export leakage
export water_loss

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
Computes evapotranspiration losses

# Arguments
- `s`: soil moisture
- `sh`: soil moisture at hygroscopic point
- `sw`: soil moisture at wilting point
- `sstar`: soil moisture below field capacity
- `emax`: maximum evapotranspiration rate
- `ew`: soil evaporation rate

"""
function evapotranspiration(s, sh, sw, sstar, emax, ew)::Float64
    if s <= sh
        et = 0.0
    elseif sh <= s <= sw
        et = ew * (s - sh) / (sw - sh)
    elseif sw <= s <= sstar
        et = ew + (emax - ew) * (s - sw) / (sstar - sw)
    else
        et = emax
    end
    return et
end


"""
Computes leakage

# Arguments
- `s`: soil moisture
- `sfc`: soil moisture at field capacity
- `b`: leakage curve exponent
- `ks`: hydralic conductivity
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
Water loss function

Combines evapotranspiration and leakage as a function of soil moisture.

# Arguments
- `s`: soil moisture
- `sh`: soil moisture at hygroscopic point
- `sw`: soil moisture at wilting point
- `sstar`: soil moisture below field capacity
- `emax`: maximum evapotranspiration rate
- `ew`: soil evaporation rate
- `sfc`: soil moisture at field capacity
- `b`: leakage curve exponent
- `ks`: hydralic conductivity
"""
function water_loss(s, sh, sw, sstar, emax, ew, sfc, b, ks)::Float64
    et = evapotranspiration(s, sh, sw, sstar, emax, ew)
    lk = leakage(s, sfc, b, ks)
    return et + lk
end

"""

"""
end
