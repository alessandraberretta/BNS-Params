import numpy as np
import bilby
import sys

# Specify the output directory and the name of the simulation.
outdir = 'outdir_bns_all'
label = 'bns_all'
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility.  This is optional!
np.random.seed(88170235)

# We are going to inject a binary neutron star waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# aligned spins of both black holes (chi_1, chi_2), etc.
injection_parameters = dict(
    mass_1=1.5, mass_2=1.3, chi_1=0.02, chi_2=0.02, luminosity_distance=360.,
    theta_jn=0.4, psi=2.659, phase=1.3, geocent_time=1735257618, 
    ra=34.9, dec=-7.1, lambda_1=400, lambda_2=450)

duration = 2000
sampling_frequency = 2 * 1024
start_time = injection_parameters['geocent_time'] + 2 - duration

# Fixed arguments passed into the source model. The analysis starts at 40 Hz.
waveform_arguments = dict(waveform_approximant='IMRPhenomPv2_NRTidal',
                          reference_frequency=50., minimum_frequency=10.0)

# Create the waveform_generator using a LAL Binary Neutron Star source function
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    waveform_arguments=waveform_arguments)

# Set up interferometers.
interferometers = bilby.gw.detector.InterferometerList(['ET'])
for interferometer in interferometers:
    interferometer.minimum_frequency = 10
interferometers.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=start_time)
interferometers.inject_signal(parameters=injection_parameters,
                              waveform_generator=waveform_generator)

# Load the default prior for binary neutron stars.
priors = bilby.gw.prior.BNSPriorDict()
priors['geocent_time'] = 1735257618
'''
priors['dec'] = -7.1
priors['ra'] =  34.9
priors['luminosity_distance'] = 360
priors['geocent_time'] = 1735257618
priors['chirp_mass'] = bilby.core.prior.Gaussian(1.215, 0.1, name='chirp_mass', unit='$M_{\\odot}$')

priors = bilby.gw.prior.BNSPriorDict()
for key in ['psi', 'geocent_time', 'ra', 'dec', 'chi_1', 'chi_2',
            'luminosity_distance', 'phase']:
    priors[key] = injection_parameters[key]
priors.pop('mass_ratio')
# priors.pop('lambda_1')
# priors.pop('lambda_2')
priors['chirp_mass'] = bilby.core.prior.Gaussian(
    1.215, 0.1, name='chirp_mass', unit='$M_{\\odot}$')
priors['symmetric_mass_ratio'] = bilby.core.prior.Uniform(
    0.1, 0.25, name='symmetric_mass_ratio')
# priors['lambda_tilde'] = bilby.core.prior.Uniform(0, 5000, name='lambda_tilde')
# priors['delta_lambda'] = bilby.core.prior.Uniform(
    # -5000, 5000, name='delta_lambda')
'''
# Initialise the likelihood by passing in the interferometer data (IFOs)
# and the waveform generator
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=interferometers, waveform_generator=waveform_generator,
    time_marginalization=False, phase_marginalization=False,
    distance_marginalization=False, priors=priors)

# Run sampler.  In this case we're going to use the `nestle` sampler
result = bilby.run_sampler(
    likelihood=likelihood, priors=priors, sampler='nestle', npoints=100,
    injection_parameters=injection_parameters, outdir=outdir, label=label,
    conversion_function=bilby.gw.conversion.generate_all_bns_parameters, 
    result_class=bilby.gw.result.CBCResult)

# result.plot_waveform_posterior(n_samples=5)

result.plot_corner()