
FUNCTION vis_duplicate, visin
  ;+
  ; Procedure:
  ;   Given an array of combined visibilities (RHESSI or STIX) with no -ve v values, vis_duplicate returns a new array with duplication
  ;   of uv points in the v<0 half-plane. The output visibility values with -ve v values are computed by converting the
  ;   corresponding visibilities with +ve v values into their conjugate values.
  ;
  ; Usage:
  ;   visout = vis_duplicate(visin)
  ;
  ; Input:
  ; visin is an array of visibility structures.
  ;
  ; Output:
  ; visout is an array of visibility structures including visibilities with -ve v values.
  ;
  ; Modification History:
  ;    Dec-10     Initial version. Anna Maria Massone
  ;    Sep-12     We preserved the Signal-to-Noise of the visibility bag. R.A. Schwartz
  ;    Sep-20     Adapted to STIX, Emma Perracchione

  
   if (visin[0].type Eq 'stx_visibility') then begin
    
    ;;;;;;;;;;;;;;;;; Check if duplication is necessary ;;;;;;;;;;;;;;;;;

    cntr_vis = N_elements(visin)

    ;;;;;;;;;;;;;;;;; If duplication is not necessary just return the input visibilities ;;;;;;;;;;;;;;;;;

    if (cntr_vis GT 60) then RETURN, visin

    ;;;;;;;;;;;;;;;;; If duplication is necessary, create a new array of visibility structures ;;;;;;;;;;;;;;;;;

    visout = replicate(stx_visibility(),2*cntr_vis)

    ;;;;;;;;;;;;;;;;; In the first half positions copy the input visibilities ;;;;;;;;;;;;;;;;;
    
    visout[0:cntr_vis-1] = visin

    ;;;;;;;;;;;;;;;;; In the second half positions duplicate the visibilities ;;;;;;;;;;;;;;;;;

    visout[cntr_vis:2*cntr_vis-1].u              = -visin.u
    visout[cntr_vis:2*cntr_vis-1].v              = -visin.v
    visout[cntr_vis:2*cntr_vis-1].obsvis         = conj(visin.obsvis)
    visout[cntr_vis:2*cntr_vis-1].isc            = visin.isc
    visout[cntr_vis:2*cntr_vis-1].ENERGY_RANGE   = visin.ENERGY_RANGE
    visout[cntr_vis:2*cntr_vis-1].time_range     = visin.time_range
    visout[cntr_vis:2*cntr_vis-1].totflux        = visin.totflux
    visout[cntr_vis:2*cntr_vis-1].sigamp         = visin.sigamp *sqrt(2.0)
    visout[0:cntr_vis-1].sigamp *= sqrt(2.0)

    visout[cntr_vis:2*cntr_vis-1].xyoffset       = visin.xyoffset
    visout[cntr_vis:2*cntr_vis-1].calibrated     = visin.calibrated
    visout[cntr_vis:2*cntr_vis-1].phase_sense    = visin.phase_sense
    visout[cntr_vis:2*cntr_vis-1].live_time      = visin.live_time
    visout[cntr_vis:2*cntr_vis-1].label      = visin.label
    endif else begin
      
      i  = hsi_vis_select(visin, PAOUT=pa)
      neg = WHERE (pa GE 180, nneg)

      ;;;;;;;;;;;;;;;;; If duplication is not necessary just return the input visibilities ;;;;;;;;;;;;;;;;;

      if (nneg GT 0) then RETURN, visin

      ;;;;;;;;;;;;;;;;; If duplication is necessary, create a new array ;;;;;;;;;;;;;;;;;

      cntr_vis = N_elements(visin)
      visout = replicate({hsi_vis},2*cntr_vis)

      ;;;;;;;;;;;;;;;;; In the first half positions copy the input visibilities ;;;;;;;;;;;;;;;;;
      
      visout[0:cntr_vis-1]=visin

      ;;;;;;;;;;;;;;;;; In the second half positions duplicate the visibilities ;;;;;;;;;;;;;;;;;

      visout[cntr_vis:2*cntr_vis-1].u          = -visin.u
      visout[cntr_vis:2*cntr_vis-1].v          = -visin.v
      visout[cntr_vis:2*cntr_vis-1].obsvis     = conj(visin.obsvis)
      visout[cntr_vis:2*cntr_vis-1].isc        = visin.isc
      visout[cntr_vis:2*cntr_vis-1].harm       = visin.harm
      visout[cntr_vis:2*cntr_vis-1].erange     = visin.erange
      visout[cntr_vis:2*cntr_vis-1].trange     = visin.trange
      visout[cntr_vis:2*cntr_vis-1].totflux    = visin.totflux
      visout[cntr_vis:2*cntr_vis-1].sigamp     = visin.sigamp *sqrt(2.0)
      visout[0:cntr_vis-1].sigamp *= sqrt(2.0)
      visout[cntr_vis:2*cntr_vis-1].chi2       = visin.chi2
      visout[cntr_vis:2*cntr_vis-1].xyoffset   = visin.xyoffset
      
    endelse
    
    return, visout
    
end