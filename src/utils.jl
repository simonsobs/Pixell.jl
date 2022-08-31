import PhysicalConstants.CODATA2018 as consts
import Unitful: ustrip

@doc raw"""
The derivative of the planck spectrum with respect to temperature, evaluated
at frequencies f and temperature T, in units of Jy/sr/K.

A blackbody has intensity ``I = 2hf^3/c^2/(\exp{hf/kT}-1) = V/(\exp{x}-1)``
with ``V = 2hf^3/c^2``, ``x = hf/kT``.

```math
dI/dx = -V/(\exp{x}-1)^2 * \exp{x}
```

```math
\begin{aligned}
dI/dT &= dI/dx * dx/dT \\
      &= 2hf^3/c^2/(\exp{x}-1)^2*\exp{x} * hf/k / T^2 \\
      &= 2*h^2*f^4/c^2/k/T^2 * \exp{x}/(\exp{x}-1)^2 \\
      &= 2*x^4 * k^3*T^2/(h^2*c^2) * \exp{x}/(\exp{x}-1)^2 \\
      &= \dots /(4* \sinh(x/2)^2)
\end{aligned}
```
"""
function dplanck(f, T=2.72548)
    c = ustrip(consts.c_0)
    k = ustrip(consts.k_B)
    h = ustrip(consts.h)

    x = h*f/(k*T)
    dIdT = 2*x^4 * k^3*T^2/(h^2*c^2) / (4*sinh(x/2)^2) * 1e26
    return dIdT
end



# utilities for applying an FFTLog transformation
"""
FFTLog from Hamilton 2000. Construct this type with `plan_fftlog`.
"""
struct FFTLogPlan{T, OT, AA<:AbstractArray{Complex{T},1},
                  AAR<:AbstractArray{T,1}, PT<:Plan, IPT<:Plan}
    L::T
    N::Int
    μ::OT
    q::T
    r₀::T
    k₀r₀::T
    uₘ::AA
    r::AAR
    k::AAR
    fftplan!::PT
    ifftplan!::IPT
end

function plan_fftlog(r::AA, μ, q=0.0, k₀r₀=1.0;
                     kropt=true) where {T, AA<:AbstractArray{T,1}}
    logrmin = log(first(r))
    logrmax = log(last(r))
    r₀ = exp((logrmin + logrmax)/2)
    @assert logrmin < logrmax
    N = length(r)
    L = logrmax - logrmin
    dlnr = L / (N - 1)
    if kropt
        k₀r₀ = k₀r₀_low_ringing(N, μ, q, L, k₀r₀)
    end
    k₀ = k₀r₀ / r₀
    Nhalf = N ÷ 2
    n = range(-Nhalf,Nhalf,length=N)
    k = reverse(k₀ .* exp.(n .* L / N))

    m = fftfreq(N, N)  # get indicies that go from [-N/2] to [N/2]
    uₘ_coeff = similar(r, Complex{T})
    for i in eachindex(m)
        uₘ_coeff[i] = uₘ(m[i], μ, q, dlnr, k₀r₀, N)
    end
    uₘ_coeff[N÷2+1] = real(uₘ_coeff[N÷2+1])  # eq 19
    fftplan! = plan_fft!(uₘ_coeff)
    ifftplan! = plan_ifft!(uₘ_coeff)
    return FFTLogPlan(L, N, μ, q, r₀, k₀r₀, uₘ_coeff, r, k, fftplan!, ifftplan!)
end


U_μ(μ, x) = exp(x * log(2.) - loggamma(0.5 * (μ + 1 - x)) + loggamma(0.5 * (μ + 1 + x)))
uₘ(m, μ, q, dlnr, k₀r₀, N) = (k₀r₀)^(-2π * im * m / (dlnr * N)) * U_μ(μ, q + 2π * im * m / (dlnr * N) )

function k₀r₀_low_ringing(N, μ, q, L, k₀r₀=1.0)
    # from pyfftlog
    dlnr = L / (N-1)
    xp = (μ + 1 + q) / 2
    xm = (μ + 1 - q) / 2
    y = π * im / 2 / dlnr
    zp = loggamma(xp + y)
    zm = loggamma(xm + y)
    arg = log(2 / k₀r₀) / dlnr + imag(zp + zm) / π
    return k₀r₀ * exp((arg - round(arg))* dlnr)
end

function mul!(Y, pl::FFTLogPlan, A)
    Y .= A
    Y .*= (pl.r).^(-pl.q)
    pl.fftplan! * Y
    Y .*= pl.uₘ
    pl.ifftplan! * Y
    Y .*= (pl.r).^(pl.q)
end

function ldiv!(Y, pl::FFTLogPlan, A)
    Y .= A
    Y .*= (pl.r).^(-pl.q)
    pl.fftplan! * Y
    Y ./= pl.uₘ
    pl.ifftplan! * Y
    Y .*= (pl.r).^(pl.q)
end


# adapted from pixell
struct RadialFourierTransform{T,A,AC,PL}
    dln::T
    l::A
    revl::A
    r::A
    pad::Int
    pl::PL
    fftbuffer::AC
end

function RadialFourierTransform(T=Float64; lrange=nothing, rrange=nothing, n=512, pad=256)

    if isnothing(lrange) && isnothing(rrange)
        lrange = (T(0.1), T(1e7))
    elseif isnothing(lrange)
        lrange = one(T) ./ reverse(rrange)
    end
    logl1, logl2 = log.(lrange)
    logl0 = (logl2 + logl1) / 2
    dlog = (logl2 - logl1) / n
    i0 = (n+1)/2 + pad
    l = @. exp(logl0 + ((1-i0):(n+2pad-i0)) * dlog)
    r = one(T) ./ reverse(l)
    fftbuffer = similar(r, complex(T))
    pl = plan_fftlog(r, 0; kropt=false)
    RadialFourierTransform(dlog, l, reverse(l), r, pad, pl, fftbuffer)
end



function real2harm(rft::RadialFourierTransform{T}, rprof::AbstractArray) where T
    fr = rprof .* rft.r
    rft_result_complex = rft.fftbuffer
    mul!(rft_result_complex, rft.pl, fr)
    rft_result_real = real.(reverse(rft_result_complex))  # allocate for the result
    return (2 * T(π)) .* rft_result_real ./ rft.l
end

function real2harm(rft::RadialFourierTransform{T}, rprof) where T
    real2harm(rft, rprof.(rft.r))
end

function harm2real(rft::RadialFourierTransform{T}, lprof::AbstractArray) where T
    fl = lprof .* rft.revl ./ (2 * T(π))
    rft_result_complex = rft.fftbuffer
    ldiv!(rft_result_complex, rft.pl, fl)
    rft_result_real = real.(rft_result_complex) # allocate for the result
    return rft_result_real ./ rft.r
end

function harm2real(rft::RadialFourierTransform{T}, lprof) where T
    harm2real(rft, lprof.(rft.revl))
end

function unpad(rft::RadialFourierTransform, arrs...)
    map(x->x[begin+rft.pad:end-rft.pad], arrs)
end
