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
; HISTORY: May 2025, Massa P., working only for detectors 3 to 10. No energy dependence is considered
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
  grid_orient_front_all = fff.o 
  grid_orient_rear_all = rrr.o
  
  grid_orient_front = fltarr(32)
  grid_orient_rear = fltarr(32)
  for subc_n=0,31 do begin
    
    ;; Front grid
    subc_n_front = fff.SC
    idx = where(subc_n_front eq (subc_n+1), count)
    if count eq 1 then begin
      
      grid_orient_front[subc_n] = grid_orient_front_all[idx]
    
    endif else begin
      
      grid_orient_front[subc_n] = grid_orient_front_all[idx[0]]
      
    endelse
    
    ;; Rear grid
    subc_n_rear = rrr.SC
    idx = where(subc_n_rear eq (subc_n+1), count)
    if count eq 1 then begin

      grid_orient_rear[subc_n] = grid_orient_rear_all[idx]

    endif else begin

      grid_orient_rear[subc_n] = grid_orient_rear_all[idx[0]]

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
  
  return, subc_transm

end
