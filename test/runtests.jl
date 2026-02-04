using Test
using SpectroscopyTools
using Unitful
using Statistics

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

    @testset "AxisType enum" begin
        @test time_axis isa AxisType
        @test wavelength_axis isa AxisType
        @test time_axis != wavelength_axis
    end

    @testset "PumpProbeData accessors" begin
        ppd = PumpProbeData(
            [1.0, 2.0, 3.0],
            ones(3, 4), ones(3, 4), zeros(3, 4),
            "20250101_120000", time_axis
        )
        @test xaxis(ppd) === ppd.time
        @test xaxis_label(ppd) == "Time (ps)"

        ppd_wl = PumpProbeData(
            [800.0, 850.0, 900.0],
            ones(3, 4), ones(3, 4), zeros(3, 4),
            "20250101_120000", wavelength_axis
        )
        @test xaxis_label(ppd_wl) == "Wavelength (nm)"
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
        result = fit_exp_decay(trace; irf_width=0.2)

        @test result isa ExpDecayIRFFit
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

        @test result isa ExpDecayIRFFit
        @test result.tau ≈ tau_true atol=1.5
        @test isnan(result.sigma)
        @test result.rsquared > 0.99
    end

    @testset "Biexponential fitting - synthetic" begin
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
        result = fit_biexp_decay(trace; irf_width=0.2)

        @test result isa BiexpDecayFit
        @test result.tau1 > 0
        @test result.tau2 > 0
        @test result.tau1 < result.tau2
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
        result2 = fit_exp_decay(trace; n_exp=2, irf_width=0.2)
        @test result2 isa MultiexpDecayFit
        @test n_exp(result2) == 2
        @test length(result2.taus) == 2
        @test length(result2.amplitudes) == 2
        @test all(result2.taus .> 0)
        @test result2.taus[1] <= result2.taus[2]
        @test result2.rsquared > 0.95

        w = weights(result2)
        @test length(w) == 2
        @test sum(w) ≈ 1.0 atol=1e-10

        # predict should work
        curve = predict(result2, trace)
        @test length(curve) == length(trace.time)
        @test all(isfinite, curve)
    end

    @testset "Global fitting - synthetic" begin
        t = collect(-2.0:0.1:30.0)
        tau_true = 5.0
        sigma_true = 0.3

        signal_esa = [SpectroscopyTools._exp_decay_irf_conv(ti, 0.8, tau_true, 0.0, sigma_true) + 0.01
                      for ti in t]
        signal_gsb = [SpectroscopyTools._exp_decay_irf_conv(ti, -0.5, tau_true, 0.0, sigma_true) - 0.005
                      for ti in t]

        trace_esa = TATrace(t, signal_esa)
        trace_gsb = TATrace(t, signal_gsb)

        result = fit_global([trace_esa, trace_gsb]; labels=["ESA", "GSB"], irf_width=0.2)

        @test result isa GlobalFitResult
        @test result.tau > 0
        @test result.tau ≈ tau_true atol=2.0
        @test length(result.amplitudes) == 2
        @test result.labels == ["ESA", "GSB"]
        @test result.rsquared > 0.95

        # predict
        curves = predict(result, [trace_esa, trace_gsb])
        @test length(curves) == 2
        @test length(curves[1]) == length(t)
    end

    @testset "predict - ExpDecayIRFFit" begin
        t = collect(-2.0:0.1:30.0)
        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, 1.0, 5.0, 0.0, 0.3) + 0.01
                  for ti in t]

        trace = TATrace(t, signal)
        result = fit_exp_decay(trace; irf_width=0.2)

        curve = predict(result, trace)
        @test length(curve) == length(t)
        @test all(isfinite, curve)

        # Also test with vector
        curve2 = predict(result, t)
        @test curve2 == curve
    end

    @testset "predict - BiexpDecayFit" begin
        t = collect(-2.0:0.1:30.0)
        signal = [SpectroscopyTools._exp_decay_irf_conv(ti, 0.5, 2.0, 0.0, 0.3) +
                  SpectroscopyTools._exp_decay_irf_conv(ti, 0.5, 15.0, 0.0, 0.3) + 0.01
                  for ti in t]

        trace = TATrace(t, signal)
        result = fit_biexp_decay(trace; irf_width=0.2)

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
        @test result.esa_center ≈ 2040.0 atol=5.0
        @test result.gsb_center ≈ 2060.0 atol=5.0
        @test result.anharmonicity > 0
        @test result.rsquared > 0.95

        y_fit = predict(result, ν)
        @test length(y_fit) == length(ν)
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
        @test time_index(times, 2.1) == 3
        @test time_index(times, 0.0) == 1
        @test time_index(times, 7.0) == 4
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

    @testset "Unit parsing - concentrations" begin
        @test parse_concentration("1.0M") == 1.0u"M"
        @test parse_concentration("500mM") == 500.0u"mM"
        @test parse_concentration("100μM") == 100.0u"μM"
        @test parse_concentration("10nM") == 10.0u"nM"
        @test parse_concentration("1pM") == 1.0u"pM"

        @test parse_concentration("1.5mol/L") == 1.5u"mol/L"
        @test parse_concentration("0.5 mol/L") == 0.5u"mol/L"

        # Keyboard substitution
        @test parse_concentration("100uM") == 100.0u"μM"

        # Whitespace
        @test parse_concentration("  1.0 M ") == 1.0u"M"
        @test parse_concentration("500 mM") == 500.0u"mM"

        # Scientific notation
        @test parse_concentration("1e-3M") == 1e-3u"M"
        @test parse_concentration("2.5e2mM") == 250.0u"mM"
        @test parse_concentration("1.5e-6M") ≈ 1.5u"μM" rtol=1e-10

        # Conversions
        c = parse_concentration("500mM")
        @test Unitful.uconvert(u"M", c) == 0.5u"M"

        @test Unitful.uconvert(u"mol/L", parse_concentration("1.0M")) == 1.0u"mol/L"

        # Errors
        @test_throws ArgumentError parse_concentration("")
        @test_throws ArgumentError parse_concentration("abc")
        @test_throws ArgumentError parse_concentration("1.0")
        @test_throws ArgumentError parse_concentration("1.0X")
        @test_throws ArgumentError parse_concentration("1.0kg")
        @test_throws ArgumentError parse_concentration("1.0s")
    end

    @testset "Unit parsing - time" begin
        @test parse_time("500fs") == 500.0u"fs"
        @test parse_time("2.5ps") == 2.5u"ps"
        @test parse_time("100ns") == 100.0u"ns"
        @test parse_time("1μs") == 1.0u"μs"
        @test parse_time("10ms") == 10.0u"ms"
        @test parse_time("1s") == 1.0u"s"

        # Keyboard substitution
        @test parse_time("1us") == 1.0u"μs"
        @test parse_time("100 us") == 100.0u"μs"

        # Whitespace
        @test parse_time("  500 fs ") == 500.0u"fs"
        @test parse_time("2.5 ps") == 2.5u"ps"

        # Scientific notation
        @test parse_time("1e3fs") == 1000.0u"fs"
        @test parse_time("1.5e-12s") ≈ 1.5u"ps" rtol=1e-10

        # Conversions
        t = parse_time("500fs")
        @test Unitful.uconvert(u"ps", t) ≈ 0.5u"ps"

        @test Unitful.uconvert(u"fs", parse_time("1ps")) == 1000.0u"fs"
        @test Unitful.uconvert(u"ns", parse_time("1000ps")) == 1.0u"ns"

        # Errors
        @test_throws ArgumentError parse_time("")
        @test_throws ArgumentError parse_time("100")
        @test_throws ArgumentError parse_time("1.0M")
        @test_throws ArgumentError parse_time("1.0nm")
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
        md = format_results(result)
        @test occursin("Peak Fit Results", md)
        @test occursin("2062.3", md)
        @test occursin("24.7", md)
        @test occursin("0.452", md)
        @test occursin("lorentzian", md)
        @test occursin("0.9985", md)
        @test occursin("2000", md)
        @test occursin("2100", md)
    end

    @testset "format_results - ExpDecayIRFFit" begin
        result = ExpDecayIRFFit(0.5, 8.3, 0.1, 0.25, 0.01, :esa, zeros(10), 0.9923)
        md = format_results(result)
        @test occursin("Exponential Decay Fit", md)
        @test occursin("ESA", md)
        @test occursin("8.3", md)
        @test occursin("0.25", md)
        @test occursin("0.9923", md)
    end

    @testset "format_results - ExpDecayFit" begin
        result = ExpDecayFit(0.5, 8.3, 0.01, 5, :gsb, zeros(10), 0.9901)
        md = format_results(result)
        @test occursin("Exponential Decay Fit", md)
        @test occursin("GSB", md)
        @test occursin("8.3", md)
        @test occursin("0.9901", md)
    end

    @testset "format_results - BiexpDecayFit" begin
        result = BiexpDecayFit(1.5, 15.0, 0.3, 0.2, 0.1, 0.25, 0.01, :esa, zeros(10), 0.9967)
        md = format_results(result)
        @test occursin("Biexponential Decay Fit", md)
        @test occursin("ESA", md)
        @test occursin("1.5", md)
        @test occursin("15.0", md)
        @test occursin("Fast", md)
        @test occursin("Slow", md)
        @test occursin("0.9967", md)
    end

    @testset "format_results - MultiexpDecayFit" begin
        result = MultiexpDecayFit(
            [0.5, 5.0, 50.0], [0.3, 0.4, 0.1],
            0.1, 0.25, 0.01, :esa, zeros(10), 0.9991
        )
        md = format_results(result)
        @test occursin("Multi-exponential Decay Fit", md)
        @test occursin("3 components", md)
        @test occursin("0.5", md)
        @test occursin("5.0", md)
        @test occursin("50.0", md)
        @test occursin("0.9991", md)
    end

    @testset "format_results - GlobalFitResult" begin
        result = GlobalFitResult(
            8.5, 0.25, 0.1,
            [0.5, -0.3], [0.01, -0.005],
            ["ESA", "GSB"],
            0.9945,
            [0.9950, 0.9940],
            [zeros(10), zeros(10)]
        )
        md = format_results(result)
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
        fwhm = irf_fwhm(sigma)
        @test fwhm ≈ 2 * sqrt(2 * log(2)) * sigma rtol=1e-10

        pfwhm = pulse_fwhm(sigma)
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

end
