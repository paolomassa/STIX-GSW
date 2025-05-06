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
; OUTPUTS:
;
;   Transmission value for every subcollimator
;
; HISTORY: May 2025, Massa P., first version (working only for detectors 3 to 10). No energy dependence is considered
;          
; CONTACT:
;   paolo.massa@fhnw.ch
;-


function stx_subc_transmission, flare_loc

  ;; Read grid orientation. It is used to compute off-axis angle
  restore,loc_file( 'grid_temp.sav', path = getenv('STX_GRID') )
  fff=read_ascii(loc_file( 'grid_param_front.txt', path = getenv('STX_GRID') ),temp=grid_temp)
  rrr=read_ascii(loc_file( 'grid_param_rear.txt', path = getenv('STX_GRID') ),temp=grid_temp)
  
  ;; Orientation of the slits of the grid as seen from the detector side
  grid_orient_front = fff.o 
  grid_orient_rear = rrr.o
  
  grid_orient_avg = (grid_orient_front + grid_orient_rear) / 2.

  ;; Off-axis angle
  flare_loc_deg = flare_loc / 3600. ;; Convert coordinates to deg
  theta = flare_loc_deg[0] * cos(grid_orient_avg * !dtor) + flare_loc_deg[1] * sin(grid_orient_avg * !dtor)

  ;; Read intercept and slope of the transmission linear fits
  fpath = loc_file( 'CFL_subcoll_transmission.txt', path = getenv('STX_GRID') )
  readcol, fpath,subc_n,subc_label,intercept,slope,$
    skipline = 1, format = 'I,A,F,F'
  
  subc_transm = intercept + slope * theta
  
  return, subc_transm

end
