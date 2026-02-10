using Test
using SpectroscopyTools
using Unitful
using Statistics
using LinearAlgebra: rank
using JSON
import Random
Random.seed!(42)

@testset "SpectroscopyTools.jl" begin

    @testset "Type hierarchy" begin
        @test TATrace <: AbstractSpectroscopyData
        @test TASpectrum <: AbstractSpectroscopyData
        @test TAMatrix <: AbstractSpectroscopyData
    end

    @testset "AbstractSpectroscopyData interface - TATrace" begin
        trace = TATrace([0.0, 1.0, 2.0], [0.1, 0.5, 0.3])

        @test xdata(trace) == [0.0, 1.0, 2.0]
        @test ydata(trace) == [0.1, 0.5, 0.3]
        @test zdata(trace) === nothing
        @test xlabel(trace) == "Time (ps)"
        @test ylabel(trace) == "ΔA"
        @test is_matrix(trace) == false
    end

    @testset "AbstractSpectroscopyData interface - TASpectrum" begin
        spec = TASpectrum([2000.0, 2050.0, 2100.0], [0.1, 0.5, 0.3])

        @test xdata(spec) == [2000.0, 2050.0, 2100.0]
        @test ydata(spec) == [0.1, 0.5, 0.3]
        @test zdata(spec) === nothing
        @test xlabel(spec) == "Wavenumber (cm⁻¹)"
        @test ylabel(spec) == "ΔA"
        @test is_matrix(spec) == false
    end

    @testset "AbstractSpectroscopyData interface - TAMatrix" begin
        time = [0.0, 1.0, 2.0]
        wavelength = [800.0, 850.0, 900.0]
        data = rand(3, 3)
        matrix = TAMatrix(time, wavelength, data)

        @test xdata(matrix) == wavelength
        @test ydata(matrix) == time
        @test zdata(matrix) === data
        @test xlabel(matrix) == "Wavelength (nm)"
        @test ylabel(matrix) == "Time (ps)"
        @test zlabel(matrix) == "ΔA"
        @test is_matrix(matrix) == true

        # Test wavenumber detection
        matrix_wn = TAMatrix(time, [1900.0, 2000.0, 2100.0], data)
        @test xlabel(matrix_wn) == "Wavenumber (cm⁻¹)"
    end

    @testset "Extended interface - source_file, npoints, title" begin
        # TATrace
        trace = TATrace([0.0, 1.0, 2.0], [0.1, 0.5, 0.3];
                        metadata=Dict{Symbol,Any}(:filename => "test.lvm"))
        @test source_file(trace) == "test.lvm"
        @test npoints(trace) == 3
        @test title(trace) == "test.lvm"

        # TATrace without filename
        trace_empty = TATrace([0.0, 1.0], [0.1, 0.2])
        @test source_file(trace_empty) == ""
        @test npoints(trace_empty) == 2
        @test title(trace_empty) == ""

        # TASpectrum
        spec = TASpectrum([2000.0, 2050.0, 2100.0], [0.1, 0.5, 0.3];
                          metadata=Dict{Symbol,Any}(:filename => "spec.lvm"))
        @test source_file(spec) == "spec.lvm"
        @test npoints(spec) == 3
        @test title(spec) == "spec.lvm"

        # TAMatrix
        time = [0.0, 1.0, 2.0]
        wavelength = [800.0, 850.0, 900.0, 950.0]
        data = rand(3, 4)
        matrix = TAMatrix(time, wavelength, data;
                          metadata=Dict{Symbol,Any}(:source => "broadband-TA/"))
        @test source_file(matrix) == "broadband-TA/"
        @test npoints(matrix) == (3, 4)
        @test title(matrix) == "broadband-TA/"
    end

    @testset "TAMatrix indexing" begin
        time = [0.0, 1.0, 2.0, 3.0, 4.0]
        wavelength = [700.0, 750.0, 800.0, 850.0]
        data = rand(5, 4)
        matrix = TAMatrix(time, wavelength, data)

        # Extract TATrace at wavelength
        trace = matrix[λ=800]
        @test trace isa TATrace
        @test length(trace.time) == 5
        @test trace.wavelength ≈ 800.0

        # Extract TASpectrum at time
        spec = matrix[t=2.0]
        @test spec isa TASpectrum
        @test length(spec.wavenumber) == 4
        @test spec.time_delay ≈ 2.0

        # Error cases
        @test_throws ErrorException matrix[]
        @test_throws ErrorException matrix[λ=800, t=1.0]
    end

    @testset "fit_peaks with raw vectors" begin
        # Synthetic single lorentzian peak
        x = collect(1900.0:0.5:2200.0)
        y = @. 0.5 / (1 + ((x - 2050.0) / 10.0)^2) + 0.01

        result = fit_peaks(x, y; n_peaks=1)
        @test result isa MultiPeakFitResult
        @test length(result) == 1
        @test result[1][:center].value ≈ 2050.0 atol=1.0
        @test result[1][:fwhm].value ≈ 20.0 atol=2.0
        @test result.r_squared > 0.99
    end

    @testset "fit_peaks multi-peak" begin
        # Synthetic two-peak spectrum
        x = collect(1900.0:0.5:2200.0)
        y = @. 0.5 / (1 + ((x - 2020.0) / 8.0)^2) + 0.3 / (1 + ((x - 2080.0) / 6.0)^2) + 0.01

        result = fit_peaks(x, y; n_peaks=2)
        @test result isa MultiPeakFitResult
        @test length(result) == 2

        centers = sort([result[1][:center].value, result[2][:center].value])
        @test centers[1] ≈ 2020.0 atol=5.0
        @test centers[2] ≈ 2080.0 atol=5.0

        y1 = predict_peak(result, 1)
        y2 = predict_peak(result, 2)
        @test length(y1) == length(x)
        @test length(y2) == length(x)

        @test result.r_squared > 0.99
    end

    @testset "MultiPeakFitResult indexing and iteration" begin
        x = collect(1900.0:0.5:2200.0)
        y = @. 0.5 / (1 + ((x - 2050.0) / 10.0)^2) + 0.01

        result = fit_peaks(x, y; n_peaks=1)
        @test result[1] isa PeakFitResult
        @test firstindex(result) == 1
        @test lastindex(result) == 1

        count = 0
        for pk in result
            @test pk isa PeakFitResult
            count += 1
        end
        @test count == length(result)
    end

    @testset "MultiPeakFitResult predict and residuals" begin
        x = collect(1900.0:0.5:2200.0)
        y = @. 0.5 / (1 + ((x - 2050.0) / 10.0)^2) + 0.01

        result = fit_peaks(x, y; n_peaks=1)

        y_fit = predict(result)
        @test length(y_fit) == result.npoints

        y_bl = predict_baseline(result)
        @test length(y_bl) == result.npoints

        res = residuals(result)
        @test length(res) == result.npoints
    end

    @testset "Baseline correction - ALS" begin
        # Synthetic spectrum: peaks on a linear baseline
        x = collect(1.0:200.0)
        baseline_true = 0.1 .+ 0.001 .* x
        peaks = 2.0 .* exp.(-((x .- 50.0) ./ 10.0).^2) .+ 1.5 .* exp.(-((x .- 150.0) ./ 8.0).^2)
        y = baseline_true .+ peaks

        bl = als_baseline(y; λ=1e5, p=0.01)
        @test length(bl) == length(y)
        # Baseline should be below the peaks
        @test maximum(bl) < maximum(y)
    end

    @testset "Baseline correction - arPLS" begin
        x = collect(1.0:200.0)
        baseline_true = 0.1 .+ 0.001 .* x
        peaks = 2.0 .* exp.(-((x .- 50.0) ./ 10.0).^2)
        y = baseline_true .+ peaks

        bl = arpls_baseline(y; λ=1e5)
        @test length(bl) == length(y)
        @test maximum(bl) < maximum(y)
    end

    @testset "Baseline correction - SNIP" begin
        x = collect(1.0:200.0)
        baseline_true = 0.5 .* ones(200)
        peaks = 3.0 .* exp.(-((x .- 100.0) ./ 15.0).^2)
        y = baseline_true .+ peaks

        bl = snip_baseline(y; iterations=40)
        @test length(bl) == length(y)
        @test maximum(bl) < maximum(y)
    end

    @testset "Baseline correction - correct_baseline API" begin
        y = randn(100) .+ 10.0
        result = correct_baseline(y; method=:arpls)
        @test haskey(result, :y)
        @test haskey(result, :baseline)
        @test length(result.y) == length(y)
        @test length(result.baseline) == length(y)

        # With x values
        x = collect(1.0:100.0)
        result2 = correct_baseline(x, y; method=:als)
        @test haskey(result2, :x)
        @test haskey(result2, :y)
        @test haskey(result2, :baseline)

        # Invalid method
        @test_throws ArgumentError correct_baseline(y; method=:invalid)
    end

    @testset "Exponential fitting - synthetic decay" begin
        # Create synthetic decay: A * exp(-(t-t0)/tau) convolved with Gaussian IRF
        t = collect(-5.0:0.1:50.0)
        tau_true = 8.0
        A_true = 1.0
        sigma_true = 0.3
        offset_true = 0.02

        # Use the internal IRF convolution to generate test data
        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, A_true, tau_true, 0.0, sigma_true) + offset_true
                  for ti in t]
        # Add small noise
        signal .+= 0.001 .* randn(length(signal))

        trace = TATrace(t, signal)
        result = fit_exp_decay(trace; irf=true, irf_width=0.2)

        @test result isa ExpDecayFit
        @test result.tau ≈ tau_true atol=1.0
        @test !isnan(result.sigma)
        @test result.rsquared > 0.99
    end

    @testset "Exponential fitting - no IRF" begin
        t = collect(0.0:0.1:50.0)
        tau_true = 10.0
        A_true = 0.8
        offset_true = 0.01

        signal = @. A_true * exp(-t / tau_true) + offset_true
        signal .+= 0.001 .* randn(length(signal))

        trace = TATrace(t, signal)
        result = fit_exp_decay(trace; irf=false)

        @test result isa ExpDecayFit
        @test result.tau ≈ tau_true atol=1.5
        @test isnan(result.sigma)
        @test result.rsquared > 0.99
    end

    @testset "Biexponential fitting via n_exp=2" begin
        t = collect(-2.0:0.1:50.0)
        tau1_true = 2.0
        tau2_true = 20.0
        A1_true = 0.4
        A2_true = 0.6
        sigma_true = 0.3

        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, A1_true, tau1_true, 0.0, sigma_true) +
                  SpectroscopyTools._exp_decay_irf_conv(ti, A2_true, tau2_true, 0.0, sigma_true) + 0.01
                  for ti in t]
        signal .+= 0.002 .* randn(length(signal))

        trace = TATrace(t, signal)
        result = fit_exp_decay(trace; n_exp=2, irf=true, irf_width=0.2)

        @test result isa MultiexpDecayFit
        @test length(result.taus) == 2
        @test all(result.taus .> 0)
        @test result.taus[1] < result.taus[2]
        @test result.rsquared > 0.95
    end

    @testset "Multi-exponential fitting (n_exp parameter)" begin
        t = collect(-2.0:0.1:50.0)
        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, 0.5, 3.0, 0.0, 0.3) +
                  SpectroscopyTools._exp_decay_irf_conv(ti, 0.3, 15.0, 0.0, 0.3) + 0.01
                  for ti in t]
        signal .+= 0.002 .* randn(length(signal))

        trace = TATrace(t, signal)

        # n_exp=2
        result2 = fit_exp_decay(trace; n_exp=2, irf=true, irf_width=0.2)
        @test result2 isa MultiexpDecayFit
        @test SpectroscopyTools.n_exp(result2) == 2
        @test length(result2.taus) == 2
        @test length(result2.amplitudes) == 2
        @test all(result2.taus .> 0)
        @test result2.taus[1] <= result2.taus[2]
        @test result2.rsquared > 0.95

        w = SpectroscopyTools.weights(result2)
        @test length(w) == 2
        @test sum(w) ≈ 1.0 atol=1e-10

        # predict should work
        curve = predict(result2, trace)
        @test length(curve) == length(trace.time)
        @test all(isfinite, curve)
    end

    @testset "Global fitting - synthetic (n_exp=1)" begin
        t = collect(-2.0:0.1:30.0)
        tau_true = 5.0
        sigma_true = 0.3

        signal_esa = [SpectroscopyTools._exp_decay_irf_conv(ti, 0.8, tau_true, 0.0, sigma_true) + 0.01
                      for ti in t]
        signal_gsb = [SpectroscopyTools._exp_decay_irf_conv(ti, -0.5, tau_true, 0.0, sigma_true) - 0.005
                      for ti in t]

        trace_esa = TATrace(t, signal_esa)
        trace_gsb = TATrace(t, signal_gsb)

        result = fit_global([trace_esa, trace_gsb]; n_exp=1, labels=["ESA", "GSB"], irf_width=0.2)

        @test result isa GlobalFitResult
        @test length(result.taus) == 1
        @test result.taus[1] > 0
        @test result.taus[1] ≈ tau_true atol=2.0
        @test size(result.amplitudes) == (2, 1)
        @test result.labels == ["ESA", "GSB"]
        @test result.rsquared > 0.95
        @test isnothing(result.wavelengths)
        @test SpectroscopyTools.n_exp(result) == 1

        # predict
        curves = predict(result, [trace_esa, trace_gsb])
        @test length(curves) == 2
        @test length(curves[1]) == length(t)
    end

    @testset "Global fitting - multi-exp (n_exp=2)" begin
        t = collect(-2.0:0.1:50.0)
        tau1_true = 2.0
        tau2_true = 15.0
        sigma_true = 0.3

        # Trace 1: both components positive (ESA-like)
        signal1 = [SpectroscopyTools._exp_decay_irf_conv(ti, 0.6, tau1_true, 0.0, sigma_true) +
                   SpectroscopyTools._exp_decay_irf_conv(ti, 0.4, tau2_true, 0.0, sigma_true) + 0.01
                   for ti in t]
        # Trace 2: both components negative (GSB-like)
        signal2 = [SpectroscopyTools._exp_decay_irf_conv(ti, -0.3, tau1_true, 0.0, sigma_true) +
                   SpectroscopyTools._exp_decay_irf_conv(ti, -0.5, tau2_true, 0.0, sigma_true) - 0.005
                   for ti in t]

        trace1 = TATrace(t, signal1)
        trace2 = TATrace(t, signal2)

        result = fit_global([trace1, trace2]; n_exp=2, irf_width=0.2, labels=["ESA", "GSB"])

        @test result isa GlobalFitResult
        @test length(result.taus) == 2
        @test result.taus[1] < result.taus[2]  # sorted fast→slow
        @test result.taus[1] ≈ tau1_true atol=2.0
        @test result.taus[2] ≈ tau2_true atol=5.0
        @test size(result.amplitudes) == (2, 2)
        @test result.rsquared > 0.95
        @test SpectroscopyTools.n_exp(result) == 2
    end

    @testset "Global fitting - TAMatrix dispatch" begin
        t = collect(-2.0:0.2:30.0)
        wavelength = collect(500.0:20.0:700.0)
        tau_true = 5.0
        sigma_true = 0.3

        data = zeros(length(t), length(wavelength))
        for (j, wl) in enumerate(wavelength)
            amp = 0.5 * sin((wl - 500) / 200 * pi)
            for (i, ti) in enumerate(t)
                data[i, j] = SpectroscopyTools._exp_decay_irf_conv(ti, amp, tau_true, 0.0, sigma_true) + 0.01
            end
        end

        matrix = TAMatrix(t, wavelength, data)
        result = fit_global(matrix; n_exp=1, irf_width=0.2)

        @test result isa GlobalFitResult
        @test !isnothing(result.wavelengths)
        @test length(result.wavelengths) == length(wavelength)
        @test result.taus[1] ≈ tau_true atol=2.0
        @test size(result.amplitudes, 1) == length(wavelength)
        @test result.rsquared > 0.95

        # DAS accessor
        d = das(result)
        @test size(d) == (1, length(wavelength))
    end

    @testset "das accessor - error without wavelengths" begin
        r = GlobalFitResult(
            [5.0], 0.25, 0.1,
            reshape([0.5, -0.3], 2, 1), [0.01, -0.005],
            ["ESA", "GSB"], nothing,
            0.9945, [0.9950, 0.9940],
            [zeros(10), zeros(10)]
        )
        @test_throws ErrorException das(r)
    end

    @testset "predict - ExpDecayFit" begin
        t = collect(-2.0:0.1:30.0)
        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, 1.0, 5.0, 0.0, 0.3) + 0.01
                  for ti in t]

        trace = TATrace(t, signal)
        result = fit_exp_decay(trace; irf=true, irf_width=0.2)

        curve = predict(result, trace)
        @test length(curve) == length(t)
        @test all(isfinite, curve)

        # Also test with vector
        curve2 = predict(result, t)
        @test curve2 == curve
    end

    @testset "predict - MultiexpDecayFit (n_exp=2)" begin
        t = collect(-2.0:0.1:30.0)
        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, 0.5, 2.0, 0.0, 0.3) +
                  SpectroscopyTools._exp_decay_irf_conv(ti, 0.5, 15.0, 0.0, 0.3) + 0.01
                  for ti in t]

        trace = TATrace(t, signal)
        result = fit_exp_decay(trace; n_exp=2, irf=true, irf_width=0.2)

        curve = predict(result, trace)
        @test length(curve) == length(t)
        @test all(isfinite, curve)
    end

    @testset "TA spectrum fitting - synthetic" begin
        ν = collect(1950.0:1.0:2150.0)
        # ESA at 2040, GSB at 2060
        esa = @. 0.005 * exp(-4 * log(2) * ((ν - 2040.0) / 15.0)^2)
        gsb = @. 0.008 * exp(-4 * log(2) * ((ν - 2060.0) / 18.0)^2)
        signal = esa .- gsb

        spec = TASpectrum(ν, signal)
        result = fit_ta_spectrum(spec; region=(1980, 2120))

        @test result isa TASpectrumFit
        @test length(result.peaks) == 2
        @test result[:esa].center ≈ 2040.0 atol=5.0
        @test result[:gsb].center ≈ 2060.0 atol=5.0
        @test result[:esa].label == :esa
        @test result[:gsb].label == :gsb
        @test anharmonicity(result) > 0
        @test result.rsquared > 0.95

        y_fit = predict(result, ν)
        @test length(y_fit) == length(ν)

        # predict_peak decomposes into individual contributions
        esa_contrib = predict_peak(result, 1, ν)
        gsb_contrib = predict_peak(result, 2, ν)
        @test all(esa_contrib .>= 0)  # ESA is positive
        @test all(gsb_contrib .<= 0)  # GSB is negative
    end

    @testset "TA spectrum fitting - three peaks" begin
        ν = collect(1900.0:1.0:2200.0)
        esa = @. 0.004 * exp(-4 * log(2) * ((ν - 2030.0) / 20.0)^2)
        gsb = @. 0.008 * exp(-4 * log(2) * ((ν - 2060.0) / 15.0)^2)
        se = @. 0.003 * exp(-4 * log(2) * ((ν - 2100.0) / 25.0)^2)
        signal = esa .- gsb .- se

        spec = TASpectrum(ν, signal)
        result = fit_ta_spectrum(spec; peaks=[:esa, :gsb, :se])

        @test length(result.peaks) == 3
        @test result[:esa].label == :esa
        @test result[:gsb].label == :gsb
        @test result[:se].label == :se
        @test result.rsquared > 0.9
    end

    @testset "TA spectrum fitting - per-peak models" begin
        ν = collect(1950.0:1.0:2150.0)
        esa = @. 0.005 * exp(-4 * log(2) * ((ν - 2040.0) / 15.0)^2)
        gsb = @. 0.008 * exp(-4 * log(2) * ((ν - 2060.0) / 18.0)^2)
        signal = esa .- gsb

        spec = TASpectrum(ν, signal)
        result = fit_ta_spectrum(spec; peaks=[(:esa, lorentzian), (:gsb, gaussian)],
                                 region=(1980, 2120))

        @test result[:esa].model == "lorentzian"
        @test result[:gsb].model == "gaussian"
        @test result.rsquared > 0.9
    end

    @testset "TA spectrum fitting - multiple peaks same type" begin
        # W(CO)6-like: 3 ESA + 3 GSB at well-separated positions
        ν = collect(1850.0:0.5:2150.0)
        esa1 = @. 0.003 * exp(-4 * log(2) * ((ν - 1960.0) / 10.0)^2)
        esa2 = @. 0.005 * exp(-4 * log(2) * ((ν - 1990.0) / 12.0)^2)
        esa3 = @. 0.002 * exp(-4 * log(2) * ((ν - 2060.0) / 8.0)^2)
        gsb1 = @. 0.004 * exp(-4 * log(2) * ((ν - 1970.0) / 10.0)^2)
        gsb2 = @. 0.007 * exp(-4 * log(2) * ((ν - 2000.0) / 12.0)^2)
        gsb3 = @. 0.003 * exp(-4 * log(2) * ((ν - 2070.0) / 8.0)^2)
        signal = (esa1 .+ esa2 .+ esa3) .- (gsb1 .+ gsb2 .+ gsb3)

        spec = TASpectrum(ν, signal)
        result = fit_ta_spectrum(spec; peaks=[:esa, :esa, :esa, :gsb, :gsb, :gsb])

        @test length(result.peaks) == 6
        @test count(p -> p.label == :esa, result.peaks) == 3
        @test count(p -> p.label == :gsb, result.peaks) == 3
        @test result.rsquared > 0.9

        # anharmonicity is NaN with multiple ESA/GSB
        @test isnan(anharmonicity(result))

        # Individual peak contributions have correct signs
        for i in 1:3
            @test all(predict_peak(result, i) .>= -1e-10)  # ESA peaks positive
        end
        for i in 4:6
            @test all(predict_peak(result, i) .<= 1e-10)   # GSB peaks negative
        end
    end

    @testset "Spectroscopy utilities - normalize" begin
        x = [1.0, 2.0, -3.0, 0.5]
        xn = normalize(x)
        @test maximum(abs.(xn)) ≈ 1.0
        @test xn[3] ≈ -1.0

        # Zero vector
        z = zeros(5)
        @test normalize(z) == zeros(5)
    end

    @testset "Spectroscopy utilities - time_index" begin
        times = [0.0, 1.0, 2.0, 5.0, 10.0]
        @test SpectroscopyTools.time_index(times, 2.1) == 3
        @test SpectroscopyTools.time_index(times, 0.0) == 1
        @test SpectroscopyTools.time_index(times, 7.0) == 4
    end

    @testset "Spectroscopy utilities - transmittance/absorbance" begin
        # Scalar
        @test transmittance_to_absorbance(0.1) ≈ 1.0
        @test transmittance_to_absorbance(50.0, percent=true) ≈ log10(2) atol=1e-6
        @test absorbance_to_transmittance(1.0) ≈ 0.1
        @test absorbance_to_transmittance(1.0, percent=true) ≈ 10.0

        # Vector
        T = [0.9, 0.5, 0.1]
        A = transmittance_to_absorbance(T)
        @test length(A) == 3
        T_back = absorbance_to_transmittance(A)
        @test T_back ≈ T atol=1e-10

        # Error
        @test_throws ArgumentError transmittance_to_absorbance(0.0)
        @test_throws ArgumentError transmittance_to_absorbance(-1.0)
    end

    @testset "Spectroscopy utilities - smooth_data" begin
        y = [0.0, 0.0, 1.0, 0.0, 0.0]
        ys = smooth_data(y; window=3)
        @test length(ys) == length(y)
        @test ys[3] < 1.0  # Smoothed peak should be lower
        @test ys[3] > 0.0
    end

    @testset "Spectroscopy utilities - calc_fwhm" begin
        x = collect(1.0:0.1:100.0)
        center = 50.0
        sigma = 5.0
        y = @. exp(-((x - center) / sigma)^2)

        result = calc_fwhm(x, y; smooth_window=1)
        @test result.peak_position ≈ center atol=0.2
        fwhm_expected = 2 * sigma * sqrt(log(2))
        @test result.fwhm ≈ fwhm_expected atol=1.0
    end

    @testset "Spectroscopy utilities - subtract_spectrum" begin
        ν = collect(1900.0:1.0:2100.0)
        y1 = @. 0.5 * exp(-((ν - 2000) / 20)^2)
        y2 = @. 0.1 * exp(-((ν - 2000) / 30)^2)

        # NamedTuple interface (raw vectors)
        result = subtract_spectrum((x=ν, y=y1), (x=ν, y=y2))
        @test result.x == ν
        @test result.y ≈ y1 .- y2

        # Typed interface (TASpectrum)
        spec1 = TASpectrum(ν, y1)
        spec2 = TASpectrum(ν, y2)
        result_typed = subtract_spectrum(spec1, spec2)
        @test result_typed.x == ν
        @test result_typed.y ≈ y1 .- y2

        # scale parameter
        result_scaled = subtract_spectrum(spec1, spec2; scale=0.5)
        @test result_scaled.y ≈ y1 .- 0.5 .* y2

        # Mismatched lengths → error with interpolate hint
        ν_short = collect(1900.0:2.0:2100.0)
        y_short = @. 0.1 * exp(-((ν_short - 2000) / 30)^2)
        @test_throws ErrorException subtract_spectrum((x=ν, y=y1), (x=ν_short, y=y_short))

        # Misaligned x-values (same length, different grids) → error with hint
        ν_shifted = ν .+ 0.5
        y_shifted = @. 0.1 * exp(-((ν_shifted - 2000) / 30)^2)
        @test_throws ErrorException subtract_spectrum((x=ν, y=y1), (x=ν_shifted, y=y_shifted))

        # Small misalignment (< 0.01) passes without error
        ν_tiny_shift = ν .+ 0.005
        y_tiny = @. 0.1 * exp(-((ν_tiny_shift - 2000) / 30)^2)
        result_close = subtract_spectrum((x=ν, y=y1), (x=ν_tiny_shift, y=y_tiny))
        @test result_close.y ≈ y1 .- y_tiny

        # interpolate=true with mismatched lengths → works
        result_interp = subtract_spectrum((x=ν, y=y1), (x=ν_short, y=y_short); interpolate=true)
        @test length(result_interp.x) == length(ν)
        @test length(result_interp.y) == length(ν)

        # interpolate=true with misaligned grids → works
        result_interp2 = subtract_spectrum((x=ν, y=y1), (x=ν_shifted, y=y_shifted); interpolate=true)
        @test result_interp2.x == ν
        @test length(result_interp2.y) == length(ν)

        # interpolate=true with scale
        result_interp_scaled = subtract_spectrum((x=ν, y=y1), (x=ν_short, y=y_short); interpolate=true, scale=0.5)
        @test length(result_interp_scaled.y) == length(ν)
    end

    @testset "Spectroscopy conversions - wavenumber/wavelength" begin
        @test wavenumber_to_wavelength(2000) ≈ 5000.0u"nm" rtol=1e-6
        @test wavenumber_to_wavelength(1000) ≈ 10000.0u"nm" rtol=1e-6
        @test wavenumber_to_wavelength(2000, output_unit=u"μm") ≈ 5.0u"μm" rtol=1e-6

        @test wavenumber_to_wavelength(2000u"cm^-1") ≈ 5000.0u"nm" rtol=1e-6

        @test wavelength_to_wavenumber(500) ≈ 20000.0u"cm^-1" rtol=1e-6
        @test wavelength_to_wavenumber(5000) ≈ 2000.0u"cm^-1" rtol=1e-6

        @test wavelength_to_wavenumber(500u"nm") ≈ 20000.0u"cm^-1" rtol=1e-6
        @test wavelength_to_wavenumber(5u"μm") ≈ 2000.0u"cm^-1" rtol=1e-6

        # Round-trip
        λ = wavenumber_to_wavelength(2000)
        @test wavelength_to_wavenumber(λ) ≈ 2000.0u"cm^-1" rtol=1e-6

        ν̃ = wavelength_to_wavenumber(800)
        @test wavenumber_to_wavelength(ν̃) ≈ 800.0u"nm" rtol=1e-6

        # Physical sanity
        @test wavelength_to_wavenumber(10.6u"μm") ≈ 943.4u"cm^-1" rtol=1e-3
        @test wavelength_to_wavenumber(400) ≈ 25000.0u"cm^-1" rtol=1e-6

        # Errors
        @test_throws ArgumentError wavenumber_to_wavelength(-100)
        @test_throws ArgumentError wavenumber_to_wavelength(0)
        @test_throws ArgumentError wavelength_to_wavenumber(-100)
        @test_throws ArgumentError wavelength_to_wavenumber(0)
    end

    @testset "Spectroscopy conversions - energy" begin
        @test wavenumber_to_energy(8065.54) ≈ 1.0u"eV" rtol=1e-4
        @test wavenumber_to_energy(2000) ≈ 0.248u"eV" rtol=1e-2
        @test wavenumber_to_energy(2000, output_unit=u"meV") ≈ 248.0u"meV" rtol=1e-2

        @test wavenumber_to_energy(2000u"cm^-1") ≈ 0.248u"eV" rtol=1e-2

        @test wavelength_to_energy(1239.84) ≈ 1.0u"eV" rtol=1e-4
        @test wavelength_to_energy(500) ≈ 2.48u"eV" rtol=1e-2

        @test wavelength_to_energy(500u"nm") ≈ 2.48u"eV" rtol=1e-2
        @test wavelength_to_energy(5u"μm") ≈ 0.248u"eV" rtol=1e-2

        @test energy_to_wavelength(1.0) ≈ 1239.84u"nm" rtol=1e-4
        @test energy_to_wavelength(2.0u"eV") ≈ 619.92u"nm" rtol=1e-4
        @test energy_to_wavelength(100u"meV") ≈ 12398.4u"nm" rtol=1e-4
        @test energy_to_wavelength(1.0, output_unit=u"μm") ≈ 1.24u"μm" rtol=1e-2

        # Round-trip
        E = wavelength_to_energy(800)
        @test energy_to_wavelength(E) ≈ 800.0u"nm" rtol=1e-6

        E2 = wavenumber_to_energy(2000)
        λ = energy_to_wavelength(E2)
        @test wavelength_to_wavenumber(λ) ≈ 2000.0u"cm^-1" rtol=1e-6

        # Physical sanity
        @test wavelength_to_energy(800) ≈ 1.55u"eV" rtol=1e-2
        @test wavelength_to_energy(532) ≈ 2.33u"eV" rtol=1e-2
        @test wavenumber_to_energy(1700, output_unit=u"meV") ≈ 211.0u"meV" rtol=1e-2

        # Errors
        @test_throws ArgumentError wavenumber_to_energy(-100)
        @test_throws ArgumentError wavenumber_to_energy(0)
        @test_throws ArgumentError wavelength_to_energy(-100)
        @test_throws ArgumentError wavelength_to_energy(0)
        @test_throws ArgumentError energy_to_wavelength(-1)
        @test_throws ArgumentError energy_to_wavelength(0)
    end

    @testset "Linewidth ↔ Decay time conversions" begin
        @test decay_time_to_linewidth(1.0) ≈ 0.658u"meV" rtol=1e-2
        @test decay_time_to_linewidth(1u"ps") ≈ 0.658u"meV" rtol=1e-2

        @test decay_time_to_linewidth(1000u"fs") ≈ 0.658u"meV" rtol=1e-2
        @test decay_time_to_linewidth(0.001u"ns") ≈ 0.658u"meV" rtol=1e-2
        @test decay_time_to_linewidth(100u"fs") ≈ 6.58u"meV" rtol=1e-2

        @test decay_time_to_linewidth(1u"ps", output_unit=u"eV") ≈ 0.000658u"eV" rtol=1e-2
        @test decay_time_to_linewidth(1u"ps", output_unit=u"cm^-1") ≈ 5.31u"cm^-1" rtol=1e-2
        @test decay_time_to_linewidth(1u"ps", output_unit=u"THz") ≈ 0.159u"THz" rtol=1e-2

        @test linewidth_to_decay_time(0.658) ≈ 1.0u"ps" rtol=1e-2
        @test linewidth_to_decay_time(0.658u"meV") ≈ 1.0u"ps" rtol=1e-2
        @test linewidth_to_decay_time(6.58u"meV") ≈ 100.0u"fs" rtol=1e-2

        @test linewidth_to_decay_time(5.31u"cm^-1") ≈ 1.0u"ps" rtol=1e-2
        @test linewidth_to_decay_time(0.159u"THz") ≈ 1.0u"ps" rtol=1e-2

        @test linewidth_to_decay_time(0.658u"meV", output_unit=u"fs") ≈ 1000.0u"fs" rtol=1e-2
        @test linewidth_to_decay_time(0.658u"meV", output_unit=u"ns") ≈ 0.001u"ns" rtol=1e-2

        # Round-trip
        τ = 2.5u"ps"
        Γ = decay_time_to_linewidth(τ)
        @test linewidth_to_decay_time(Γ) ≈ τ rtol=1e-10

        Γ2 = 10u"meV"
        τ2 = linewidth_to_decay_time(Γ2)
        @test decay_time_to_linewidth(τ2) ≈ Γ2 rtol=1e-10

        # Physical sanity
        @test decay_time_to_linewidth(100u"fs") ≈ 6.58u"meV" rtol=1e-2
        @test decay_time_to_linewidth(1u"ns") ≈ 0.658u"μeV" rtol=1e-2
        @test decay_time_to_linewidth(10u"fs") ≈ 65.8u"meV" rtol=1e-2

        # Errors
        @test_throws ArgumentError decay_time_to_linewidth(-1.0)
        @test_throws ArgumentError decay_time_to_linewidth(0.0)
        @test_throws ArgumentError linewidth_to_decay_time(-1.0)
        @test_throws ArgumentError linewidth_to_decay_time(0.0)
    end

    @testset "format_results - MultiPeakFitResult" begin
        pk = PeakFitResult(
            [:amplitude, :center, :fwhm],
            [0.452, 2062.3, 24.7],
            [0.002, 0.12, 0.31],
            [(0.448, 0.456), (2062.06, 2062.54), (24.08, 25.32)],
            0.9985, 0.001, 0.0001, 100,
            (2000.0, 2100.0), "lorentzian", "NH4SCN_DMF_1M"
        )
        result = MultiPeakFitResult(
            [pk], [0.01], 0, 0.9985, 0.001, 0.0001, 100,
            (2000.0, 2100.0), "NH4SCN_DMF_1M",
            [0.452, 2062.3, 24.7, 0.01], lorentzian, 3,
            collect(2000.0:1.0:2099.0), zeros(100)
        )
        md = SpectroscopyTools.format_results(result)
        @test occursin("Peak Fit Results", md)
        @test occursin("2062.3", md)
        @test occursin("24.7", md)
        @test occursin("0.452", md)
        @test occursin("lorentzian", md)
        @test occursin("0.9985", md)
        @test occursin("2000", md)
        @test occursin("2100", md)
    end

    @testset "format_results - ExpDecayFit" begin
        result = ExpDecayFit(0.5, 8.3, 0.1, 0.25, 0.01, :esa, zeros(10), 0.9923)
        md = SpectroscopyTools.format_results(result)
        @test occursin("Exponential Decay Fit", md)
        @test occursin("ESA", md)
        @test occursin("8.3", md)
        @test occursin("0.25", md)
        @test occursin("0.9923", md)
    end

    @testset "format_results - MultiexpDecayFit" begin
        result = MultiexpDecayFit(
            [0.5, 5.0, 50.0], [0.3, 0.4, 0.1],
            0.1, 0.25, 0.01, :esa, zeros(10), 0.9991
        )
        md = SpectroscopyTools.format_results(result)
        @test occursin("Multi-exponential Decay Fit", md)
        @test occursin("3 components", md)
        @test occursin("0.5", md)
        @test occursin("5.0", md)
        @test occursin("50.0", md)
        @test occursin("0.9991", md)
    end

    @testset "format_results - GlobalFitResult" begin
        result = GlobalFitResult(
            [8.5], 0.25, 0.1,
            reshape([0.5, -0.3], 2, 1), [0.01, -0.005],
            ["ESA", "GSB"], nothing,
            0.9945,
            [0.9950, 0.9940],
            [zeros(10), zeros(10)]
        )
        md = SpectroscopyTools.format_results(result)
        @test occursin("Global Fit Results", md)
        @test occursin("2 traces", md)
        @test occursin("8.5", md)
        @test occursin("ESA", md)
        @test occursin("GSB", md)
        @test occursin("0.9945", md)
        @test occursin("Shared Parameters", md)
        @test occursin("Per-Trace", md)
    end

    @testset "irf_fwhm and pulse_fwhm" begin
        sigma = 0.3
        fwhm = SpectroscopyTools.irf_fwhm(sigma)
        @test fwhm ≈ 2 * sqrt(2 * log(2)) * sigma rtol=1e-10

        pfwhm = SpectroscopyTools.pulse_fwhm(sigma)
        @test pfwhm ≈ fwhm / sqrt(2) rtol=1e-10
    end

    @testset "PeakInfo struct" begin
        pi = PeakInfo(2050.0, 0.5, 0.4, 20.0, (2040.0, 2060.0), 100)
        @test pi.position == 2050.0
        @test pi.intensity == 0.5
        @test pi.prominence == 0.4
        @test pi.width == 20.0
        @test pi.bounds == (2040.0, 2060.0)
        @test pi.index == 100

        # Show methods should not error
        io = IOBuffer()
        show(io, pi)
        @test length(String(take!(io))) > 0
    end

    @testset "peak_table" begin
        peaks = [
            PeakInfo(2050.0, 0.5, 0.4, 20.0, (2040.0, 2060.0), 100),
            PeakInfo(2080.0, 0.3, 0.2, 15.0, (2072.5, 2087.5), 160)
        ]
        table = peak_table(peaks)
        @test occursin("Position", table)
        @test occursin("2050.0", table)
        @test occursin("2080.0", table)

        @test peak_table(PeakInfo[]) == "No peaks detected"
    end

    @testset "Re-exports from CurveFit" begin
        # These should be accessible
        @test isdefined(SpectroscopyTools, :solve)
        @test isdefined(SpectroscopyTools, :NonlinearCurveFitProblem)
        @test isdefined(SpectroscopyTools, :coef)
        @test isdefined(SpectroscopyTools, :stderror)
    end

    @testset "Re-exports from CurveFitModels" begin
        @test isdefined(SpectroscopyTools, :gaussian)
        @test isdefined(SpectroscopyTools, :lorentzian)
        @test isdefined(SpectroscopyTools, :single_exponential)
        @test isdefined(SpectroscopyTools, :pseudo_voigt)
    end

    @testset "Chirp correction" begin

        @testset "subtract_background" begin
            n_time = 50
            n_wl = 20
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))
            background = repeat([0.01 * j for j in 1:n_wl]', n_time)

            signal = zeros(n_time, n_wl)
            for i in eachindex(time)
                if time[i] > 0
                    for j in eachindex(wavelength)
                        signal[i, j] = 0.1 * exp(-time[i] / 3.0) * sin(j * 0.5)
                    end
                end
            end

            data = signal .+ background
            metadata = Dict{Symbol,Any}(:source => "test")
            matrix = TAMatrix(time, wavelength, data, metadata)

            corrected = subtract_background(matrix)
            @test corrected isa TAMatrix
            @test size(corrected.data) == size(matrix.data)
            @test corrected.metadata[:background_subtracted] == true
            @test haskey(corrected.metadata, :baseline_t_range)

            pre_pump_mask = corrected.time .< -1.0
            pre_pump_data = corrected.data[pre_pump_mask, :]
            @test maximum(abs.(pre_pump_data)) < 0.01

            @test matrix.data !== corrected.data

            corrected2 = subtract_background(matrix; t_range=(-5.0, -2.0))
            @test corrected2.metadata[:baseline_t_range] == (-5.0, -2.0)
        end

        @testset "detect_chirp on synthetic data" begin
            n_time = 200
            n_wl = 100
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            ref_λ = 600.0
            chirp_fn(λ) = 0.0001 * (λ - ref_λ)^2 - 0.002 * (λ - ref_λ)

            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                λ = wavelength[j]
                t_onset = chirp_fn(λ)
                for i in eachindex(time)
                    t = time[i]
                    if t > t_onset
                        data[i, j] = 0.5 * exp(-(t - t_onset) / 3.0)
                    end
                end
            end

            metadata = Dict{Symbol,Any}(:source => "synthetic")
            matrix = TAMatrix(time, wavelength, data, metadata)

            cal = detect_chirp(matrix; order=2, reference=ref_λ, smooth_window=7, bin_width=4)
            @test cal isa ChirpCalibration
            @test cal.poly_order == 2
            @test cal.reference_λ == ref_λ
            @test cal.r_squared > 0.9
            @test length(cal.wavelength) > 0
            @test length(cal.time_offset) == length(cal.wavelength)

            poly = polynomial(cal)
            @test abs(poly(ref_λ)) < 0.5

            @test report(cal) === nothing
        end

        @testset "correct_chirp" begin
            n_time = 100
            n_wl = 50
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            ref_λ = 600.0
            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                t_onset = 0.02 * (wavelength[j] - ref_λ)
                for i in eachindex(time)
                    if time[i] > t_onset
                        data[i, j] = exp(-(time[i] - t_onset) / 2.0)
                    end
                end
            end

            metadata = Dict{Symbol,Any}(:source => "synthetic")
            matrix = TAMatrix(time, wavelength, data, metadata)

            cal = ChirpCalibration(
                collect(wavelength),
                [0.02 * (λ - ref_λ) for λ in wavelength],
                [-0.02 * ref_λ, 0.02],
                1,
                ref_λ,
                1.0,
                Dict{Symbol,Any}()
            )

            corrected = correct_chirp(matrix, cal)
            @test corrected isa TAMatrix
            @test size(corrected.data) == size(matrix.data)
            @test corrected.metadata[:chirp_corrected] == true

            inner = n_wl ÷ 4 : 3 * n_wl ÷ 4
            peak_times = [time[argmax(corrected.data[:, j])] for j in inner]
            @test std(peak_times) < 1.0
        end

        @testset "save_chirp and load_chirp round-trip" begin
            cal = ChirpCalibration(
                [500.0, 600.0, 700.0],
                [-1.0, 0.0, 1.5],
                [0.1, -0.002, 0.00001],
                2,
                600.0,
                0.998,
                Dict{Symbol,Any}(:order => 2, :smooth_window => 15)
            )

            tmpfile = tempname() * ".json"
            save_chirp(tmpfile, cal)
            @test isfile(tmpfile)

            cal2 = load_chirp(tmpfile)
            @test cal2 isa ChirpCalibration
            @test cal2.wavelength ≈ cal.wavelength
            @test cal2.time_offset ≈ cal.time_offset
            @test cal2.poly_coeffs ≈ cal.poly_coeffs
            @test cal2.poly_order == cal.poly_order
            @test cal2.reference_λ ≈ cal.reference_λ
            @test cal2.r_squared ≈ cal.r_squared
            @test cal2.metadata[:order] == 2

            rm(tmpfile)
        end

        @testset "ChirpCalibration polynomial" begin
            cal = ChirpCalibration(
                Float64[], Float64[],
                [1.0, 0.5, 0.01],
                2, 0.0, 1.0,
                Dict{Symbol,Any}()
            )

            poly = polynomial(cal)
            @test poly(0.0) ≈ 1.0
            @test poly(1.0) ≈ 1.51
            @test poly(10.0) ≈ 1.0 + 5.0 + 1.0
        end

        @testset "Chirp exports available" begin
            @test isdefined(SpectroscopyTools, :ChirpCalibration)
            @test isdefined(SpectroscopyTools, :detect_chirp)
            @test isdefined(SpectroscopyTools, :correct_chirp)
            @test isdefined(SpectroscopyTools, :subtract_background)
            @test isdefined(SpectroscopyTools, :save_chirp)
            @test isdefined(SpectroscopyTools, :load_chirp)
            @test isdefined(SpectroscopyTools, :polynomial)
        end

        @testset "subtract_background with flat matrix" begin
            n_time = 20
            n_wl = 10
            time = collect(range(-2.0, 5.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))
            data = zeros(n_time, n_wl)
            metadata = Dict{Symbol,Any}(:source => "flat")
            matrix = TAMatrix(time, wavelength, data, metadata)

            corrected = subtract_background(matrix)
            @test all(corrected.data .== 0.0)
        end

        @testset "detect_chirp :threshold method" begin
            n_time = 200
            n_wl = 80
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            ref_λ = 600.0
            chirp_fn(λ) = 0.005 * (λ - ref_λ)

            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                t_onset = chirp_fn(wavelength[j])
                for i in eachindex(time)
                    if time[i] > t_onset
                        data[i, j] = 0.5 * exp(-(time[i] - t_onset) / 3.0)
                    end
                end
            end

            matrix = TAMatrix(time, wavelength, data)

            cal = detect_chirp(matrix; method=:threshold, order=1,
                               reference=ref_λ, smooth_window=7, bin_width=4)
            @test cal isa ChirpCalibration
            @test cal.r_squared > 0.8
            @test length(cal.wavelength) > 0
        end

        @testset "detect_chirp input validation" begin
            n_time = 50
            n_wl = 20
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))
            data = rand(n_time, n_wl)
            matrix = TAMatrix(time, wavelength, data)

            @test_throws ArgumentError detect_chirp(matrix; order=0)
            @test_throws ArgumentError detect_chirp(matrix; bin_width=0)
            @test_throws ArgumentError detect_chirp(matrix; bin_width=n_wl + 1)
            @test_throws ArgumentError detect_chirp(matrix; min_signal=0.0)
            @test_throws ArgumentError detect_chirp(matrix; min_signal=1.5)
            @test_throws ArgumentError detect_chirp(matrix; threshold=-1.0)
            @test_throws ArgumentError detect_chirp(matrix; method=:invalid)
            @test_throws ArgumentError detect_chirp(matrix; method=:threshold, onset_frac=0.0)
            @test_throws ArgumentError detect_chirp(matrix; method=:threshold, onset_frac=1.0)
        end

        @testset "detect_chirp even smooth_window accepted" begin
            n_time = 200
            n_wl = 40
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                for i in eachindex(time)
                    if time[i] > 0
                        data[i, j] = 0.5 * exp(-time[i] / 3.0)
                    end
                end
            end

            matrix = TAMatrix(time, wavelength, data)
            cal = detect_chirp(matrix; smooth_window=10, order=1, bin_width=4)
            @test cal isa ChirpCalibration
            @test cal.metadata[:smooth_window] == 11
        end

        @testset "detect_chirp recovers known linear chirp" begin
            n_time = 200
            n_wl = 80
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            ref_λ = 600.0
            slope = 0.01
            chirp_fn(λ) = slope * (λ - ref_λ)

            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                t_onset = chirp_fn(wavelength[j])
                for i in eachindex(time)
                    if time[i] > t_onset
                        data[i, j] = 0.5 * exp(-(time[i] - t_onset) / 3.0)
                    end
                end
            end

            matrix = TAMatrix(time, wavelength, data)
            cal = detect_chirp(matrix; order=1, reference=ref_λ, smooth_window=7, bin_width=4)

            poly = polynomial(cal)
            @test abs(poly(ref_λ)) < 0.5
            @test abs(poly(ref_λ + 50) - slope * 50) < 1.0
            @test abs(poly(ref_λ - 50) + slope * 50) < 1.0
        end

        @testset "correct_chirp tighter tolerance" begin
            n_time = 200
            n_wl = 50
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            ref_λ = 600.0
            slope = 0.02
            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                t_onset = slope * (wavelength[j] - ref_λ)
                for i in eachindex(time)
                    if time[i] > t_onset
                        data[i, j] = exp(-(time[i] - t_onset) / 2.0)
                    end
                end
            end

            matrix = TAMatrix(time, wavelength, data)

            # Exact calibration — no detection noise
            cal = ChirpCalibration(
                collect(wavelength),
                [slope * (λ - ref_λ) for λ in wavelength],
                [-slope * ref_λ, slope],
                1, ref_λ, 1.0, Dict{Symbol,Any}()
            )

            corrected = correct_chirp(matrix, cal)
            inner = n_wl ÷ 4 : 3 * n_wl ÷ 4
            peak_times = [time[argmax(corrected.data[:, j])] for j in inner]
            @test std(peak_times) < 0.3
        end

        @testset "ChirpCalibration show methods" begin
            cal = ChirpCalibration(
                [500.0, 600.0, 700.0], [-1.0, 0.0, 1.5],
                [0.1, -0.002, 0.00001], 2, 600.0, 0.998,
                Dict{Symbol,Any}(:order => 2)
            )

            # Compact show
            io = IOBuffer()
            show(io, cal)
            compact = String(take!(io))
            @test occursin("ChirpCalibration", compact)
            @test occursin("order 2", compact)
            @test occursin("3 points", compact)

            # MIME show
            io = IOBuffer()
            show(io, MIME("text/plain"), cal)
            full = String(take!(io))
            @test occursin("Polynomial order", full)
            @test occursin("R²", full)
            @test occursin("Coefficients", full)
        end

        @testset "load_chirp with malformed JSON" begin
            tmpfile = tempname() * ".json"

            # Missing required keys
            open(tmpfile, "w") do io
                JSON.print(io, Dict("wavelength" => [1.0], "poly_coeffs" => [0.1]))
            end
            @test_throws ArgumentError load_chirp(tmpfile)
            rm(tmpfile)
        end

        @testset "metadata key rename :mad_threshold" begin
            n_time = 200
            n_wl = 80
            time = collect(range(-5.0, 15.0, length=n_time))
            wavelength = collect(range(500.0, 700.0, length=n_wl))

            data = zeros(n_time, n_wl)
            for j in eachindex(wavelength)
                for i in eachindex(time)
                    if time[i] > 0
                        data[i, j] = 0.5 * exp(-time[i] / 3.0)
                    end
                end
            end

            matrix = TAMatrix(time, wavelength, data)
            cal = detect_chirp(matrix; order=1, bin_width=4)

            @test haskey(cal.metadata, :mad_threshold)
            @test !haskey(cal.metadata, :threshold)
        end

    end

    @testset "SVD filtering" begin
        time = collect(range(-1.0, 10.0, length=50))
        wavelength = collect(range(400.0, 700.0, length=30))

        # Build a rank-2 signal: two spectral components with different kinetics
        signal = zeros(50, 30)
        for j in eachindex(wavelength)
            for i in eachindex(time)
                t = time[i]
                if t > 0
                    # Component 1: fast decay, peaks at 500 nm
                    signal[i, j] += 0.5 * exp(-t / 1.0) * exp(-((wavelength[j] - 500) / 30)^2)
                    # Component 2: slow decay, peaks at 600 nm
                    signal[i, j] += 0.3 * exp(-t / 5.0) * exp(-((wavelength[j] - 600) / 40)^2)
                end
            end
        end

        noise = 0.02 * randn(50, 30)
        noisy_data = signal .+ noise

        @testset "singular_values - TAMatrix" begin
            matrix = TAMatrix(time, wavelength, noisy_data)
            sv = singular_values(matrix)
            @test length(sv) == min(50, 30)
            @test issorted(sv, rev=true)
            @test sv[1] > sv[end]
        end

        @testset "singular_values - raw matrix" begin
            sv = singular_values(noisy_data)
            @test length(sv) == min(50, 30)
            @test issorted(sv, rev=true)
        end

        @testset "svd_filter - TAMatrix denoising" begin
            matrix = TAMatrix(time, wavelength, noisy_data)
            filtered = svd_filter(matrix; n_components=2)

            @test filtered isa TAMatrix
            @test size(filtered.data) == size(matrix.data)
            @test filtered.time == matrix.time
            @test filtered.wavelength == matrix.wavelength
            @test filtered.metadata[:svd_filtered] == true
            @test filtered.metadata[:svd_n_components] == 2

            # Filtered should be closer to true signal than noisy data
            error_noisy = sum((noisy_data .- signal).^2)
            error_filtered = sum((filtered.data .- signal).^2)
            @test error_filtered < error_noisy
        end

        @testset "svd_filter - raw matrix" begin
            filtered = svd_filter(time, wavelength, noisy_data; n_components=2)
            @test size(filtered) == size(noisy_data)

            error_noisy = sum((noisy_data .- signal).^2)
            error_filtered = sum((filtered .- signal).^2)
            @test error_filtered < error_noisy
        end

        @testset "svd_filter - n_components=1 keeps dominant component" begin
            matrix = TAMatrix(time, wavelength, noisy_data)
            filtered = svd_filter(matrix; n_components=1)

            # Rank-1 approximation
            @test rank(filtered.data) == 1
        end

        @testset "svd_filter - full rank preserves data" begin
            matrix = TAMatrix(time, wavelength, noisy_data)
            max_comp = min(size(noisy_data)...)
            filtered = svd_filter(matrix; n_components=max_comp)

            @test filtered.data ≈ noisy_data atol=1e-10
        end

        @testset "svd_filter - input validation" begin
            matrix = TAMatrix(time, wavelength, noisy_data)

            @test_throws ArgumentError svd_filter(matrix; n_components=0)
            @test_throws ArgumentError svd_filter(matrix; n_components=100)
            @test_throws DimensionMismatch svd_filter(time, wavelength[1:5], noisy_data; n_components=2)
        end
    end

end
