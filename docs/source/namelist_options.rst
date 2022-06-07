Stochastic Physics Namelist 
===========================

General options 
"""""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "NEW_LSCALE", "Recommended, set to true if use the correct calculation for decorrelation length scale."
   "NTRUNC", "Optional, Spectral resolution (e.g. T126) of random patterns, default is for model to determine proper truncation"
   "LAT_S", "Optional, number of latitude points for the gaussian grid  (must be even), default is for model to determine gaussian grid"
   "LON_S", "Optional, number of longitude points for the gaussian grid (recommend 2xLAT_S, default is for model to determine gaussian grid"
   "STOCHINI", "Optional, set to true if wanting to read in a previous random pattern"

SPPT options 
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "SPPT", "Amplitudes of random patterns."
   "SPPT_TAU", "Decorrelation timescales in seconds."
   "SPPT_LSCALE", "Decorrelation spatial scales in meters."
   "SPPT_LOGIT", "Should be true to limit the SPPT perturbations between 0 and 2.  Otherwise model will crash."
   "ISEED_SPPT", "Seeds for setting the random number sequence (ignored if stochini is true)"
   "SPPT_SIGTOP1", "lower sigma level to taper perturbations to zero (default is 0.1)"
   "SPPT_SIGTOP2", "upper sigma level to taper perturbations to zero (default is 0.025)"
   "SPPT_SFCLIMIT", ".T.=tapers the SPPT perturbations to zero at model’s lowest level (helps reduce model crashes)"
   "SPPTINT", "Optional, interval in seconds to update random pattern.  Perturbations still get applied every time-step"
   "USE_ZMTNBLCK", ".T.=do not apply perturbations below the dividing streamline that is diagnosed by the gravity wave drag, mountain blocking scheme"
   "PERT_MP", ".T.=apply SPPT perturbations to all microphysics specis. If .F. the SPPT is only applied to u,v,t,qv"
   "PERT_RADTEND", ".T.=apply SPPT perturbations to cloudy sky radiation tendencies. If .F. then do not perturb any radiative tendencies"
   "PERT_CLDS", ".T.=apply SPPT perturbations to fraction (only works for RRTMG radiation),  if using this option set PERT_RADTEND=.F."


SHUM options 
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "SHUM", "Amplitudes of random pattern."
   "SHUM_TAU", "Decorrelation timescales in seconds."
   "SHUM_LSCALE", "ecorrelation spatial scales in meters."
   "SHUM_SIGEFOLD", "e-folding lengthscale (in units of sigma) of specific humidity perturbations, default is 0.2)"
   "SHUMINT", "Optional, interval in seconds to update random pattern.  Perturbations still get applied every time-step"
   "ISEED_SHUM", "Seeds for setting the random number sequence (ignored if stochini is true)."

SKEB options
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "SKEB", "Amplitudes of random patterns."
   "SKEB_TAU", "Decorrelation timescales in seconds"
   "SKEB_LSCALE", "Decorrelation spatial scales in meters  (250)"
   "ISEED_SKEB", "Seeds for setting the random number sequence (ignored if stochini is true)."
   "SKEBNORM", "0-random pattern is stream function, 1-pattern is K.E. norm, 2-pattern is vorticity (default is 0)"
   "SKEB_VARSPECT_OPT", "0-gaussian (default), 1-power law (not tested)"
   "SKEB_NPASS", "number of passes of the del2 smoothing for the dissipation estimate (default is 11, minimum is 3)"
   "SKEB_VDOF", "the number of degrees of freedom in the vertical for the SKEB random pattern (default is 5)"
   "SKEB_SIGTOP1", "lower sigma level to taper perturbations to zero (default is 0.1)"
   "SKEB_SIGTOP2", "upper sigma level to taper perturbations to zero (0.025)"
   "SKEBINT", "Optional, interval in seconds to update random pattern.  Perturbations still get applied every time-step"

SPP options
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "SPP_VAR_LIST", "List of parameterizations to be perturbed. Check compns_stochy.F90 for options."
   "SPP_PRT_LIST", "SPP perturbation magnitudes used in each parameterization."
   "SPP_TAU", "Decorrelation timescales in seconds."
   "SPP_LSCALE", "Decorrelation spatial scales in meters."
   "ISEED_SPP", "Seeds for setting the random number sequence (ignored if stochini is true)."
   "SPP_SIGTOP1", "Lower sigma level to taper perturbations to zero (default is 0.1)"
   "SPP_SIGTOP2", "Upper sigma level to taper perturbations to zero (0.025)"
   "SPP_STDDEV_CUTOFF", "Limit for possible perturbation values in standard deviations from the mean."


Land perturbation options
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "LNDP_TYPE", "0, no perturbations. 1, old scheme (Gehne et al. 2019); 2, new scheme (Draper 2021)"
   "LNDP_VAR_LIST", "List of land perturbations parameters. Check compns_stochy.F90 for options"
   "LNDP_PRT_LIST", "Perturbation magnitudes used for each parameter perturbations."
   "LNDP_TAU", "Decorrelation timescales in seconds."
   "LNDP_LSCALE", "Decorrelation spatial scales in meters."
   "ISEED_LNDP", "Seeds for setting the random number sequence (ignored if stochini is true)."

