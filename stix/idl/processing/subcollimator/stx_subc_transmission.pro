;+
;
; NAME:
;
;   stx_subc_transmission
;
; PURPOSE:
;
;   Compute the transmission of a STIX subcollimator using a linear fit of CFL-calibrated data
;
; CALLING SEQUENCE:
;
;   subc_transmission = stx_subc_transmission(flare_loc)
;
; INPUTS:
;
;   flare_loc: bidimensional array containing the X and Y coordinate of the flare location
;             (arcsec, in the STIX coordinate frame)
;             
;   ph_in: an array of photon energies [keV] (restricted to range 1-1000)
;
; OUTPUTS:
;
;   Transmission value for every subcollimator
;
; KEYWORDS:
;   
;   ds: float. Width of the wedge model (in mm)
;   
;   dh: float. Hight of the wedge model (in mm)
;   
;   simple_transm: if set, the energy-dependence of the grid transmission is not considered
;   
;   silent: if set, no message is displayed
;
; HISTORY: May 2025, Massa P., based on the previous version by ECMD.
;                    Working only for detectors 3 to 10. 
;          
; CONTACT:
;   paolo.massa@fhnw.ch
;-


function stx_subc_transmission, flare_loc, ph_in, ds=ds, dh=dh, simple_transm = simple_transm, silent = silent

  ;; Parameters of the wedge grid transmission model. The values are derived from fitting of STIX observation
  default, ds, 5e-3
  default, dh, 5e-2

  ; To determine the transmission through the tungsten slats a Linear Attenuation Coefficient [mm-1]
  ; (1/absorption length in mm) is estimated for an each expected incoming photon energy and passed to stx_grid_transmission.
  ; The mass attenuation coefficient [cm^2/gm] is calculated using the xcom tabulated values in xsec.pro
  ; and a value of 19.3 g/cm3 is used for the density of tungsten. This is then divided by 10 to convert from [cm-1] to [mm-1].
  ;
  ; For backwards compatibility if no photon energy array is passed in the transmission is calculated for 1 keV
  ; (the lowest tabulated value). This should provide a reasonable approximation to the previously assumed
  ; fully opaque grids.

  if ~keyword_set(ph_in) then begin
    ph_in = 1.
    if ~keyword_set(silent) then begin
      if ~keyword_set(simple_transm) then message, 'No photon energies passed, calculating low energy approximation at 1 keV.', /info $
      else message, 'Simple grid transmission selected, calculating opaque approximation.', /info
    endif
  endif


  ;; Read grid orientation. It is used to compute off-axis angle
  restore,loc_file( 'grid_temp.sav', path = getenv('STX_GRID') )
  fff=read_ascii(loc_file( 'grid_param_front.txt', path = getenv('STX_GRID') ),temp=grid_temp)
  rrr=read_ascii(loc_file( 'grid_param_rear.txt', path = getenv('STX_GRID') ),temp=grid_temp)
  
  ;; Orientation of the slits of the grid as seen from the detector side
  grid_orient_front_all = fff.o 
  grid_orient_rear_all = rrr.o
  
  pitch_front_all = fff.P
  pitch_rear_all = rrr.P
  thickness_front_all = fff.THICK
  thickness_rear_all = rrr.THICK

  subc_n_front = fff.SC
  subc_n_rear = rrr.SC
  
  grid_orient_front = fltarr(32)
  grid_orient_rear = fltarr(32)
  pitch_front = fltarr(32)
  pitch_rear = fltarr(32)
  thickness_front = fltarr(32)
  thickness_rear = fltarr(32)
  
  for subc_n=0,31 do begin
    
    ;; Front grid
    subc_n_front = fff.SC
    idx = where(subc_n_front eq (subc_n+1), count)
    if count eq 1 then begin
      
      grid_orient_front[subc_n] = grid_orient_front_all[idx]
      pitch_front[subc_n] = pitch_front_all[idx]
      thickness_front[subc_n] = thickness_front_all[idx]
    
    endif else begin
      
      grid_orient_front[subc_n] = grid_orient_front_all[idx[0]]
      pitch_front[subc_n] = pitch_front_all[idx[0]]
      thickness_front[subc_n] = thickness_front_all[idx[0]]
      
    endelse
    
    ;; Rear grid
    subc_n_rear = rrr.SC
    idx = where(subc_n_rear eq (subc_n+1), count)
    if count eq 1 then begin

      grid_orient_rear[subc_n] = grid_orient_rear_all[idx]
      pitch_rear[subc_n] = pitch_rear_all[idx]
      thickness_rear[subc_n] = thickness_rear_all[idx]

    endif else begin

      grid_orient_rear[subc_n] = grid_orient_rear_all[idx[0]]
      pitch_rear[subc_n] = pitch_rear_all[idx[0]]
      thickness_rear[subc_n] = thickness_rear_all[idx[0]]

    endelse
    
  endfor  
  
  grid_orient_avg = (grid_orient_front + grid_orient_rear) / 2.

  ;; Off-axis angle
  flare_loc_deg = flare_loc / 3600. ;; Convert coordinates to deg
  theta = flare_loc_deg[0] * cos(grid_orient_avg * !dtor) + flare_loc_deg[1] * sin(grid_orient_avg * !dtor)

  ;; Read intercept and slope of the transmission linear fits
  fpath = loc_file( 'CFL_subcoll_transmission.txt', path = getenv('STX_GRID') )
  readcol, fpath,subc_n,subc_label,intercept,slope,$
    skipline = 1, format = 'I,A,F,F', /silent
  
  subc_transm = intercept + slope * theta
  
  if ~keyword_set(simple_transm) then begin
    
    ;; Compute slit-to-pitch ratio as the sqrt of the subcoll. transmission
    slit_to_pitch = sqrt(subc_transm)

    slit_front = slit_to_pitch*pitch_front
    slit_rear = slit_to_pitch*pitch_rear

    ;; Compute absorption
    mass_attenuation = xsec(ph_in, (Element2Z('W'))[0], 'AB', /cm2perg, error=error, use_xcom=1)
    gmcm = 19.30
    L = 1. / (mass_attenuation*gmcm/10.)

    transm_front = stx_grid_transmission(pitch_front, slit_front, thickness_front, L, ds, dh)
    transm_rear = stx_grid_transmission(pitch_rear, slit_rear, thickness_rear, L, ds, dh)
    
    return, transm_front * transm_rear
  
  endif else begin
    
    return, subc_transm
    
  endelse

end
