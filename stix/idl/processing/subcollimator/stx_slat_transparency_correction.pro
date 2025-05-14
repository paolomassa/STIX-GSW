;+
;
; NAME:
;
;   stx_slat_transparency_correction
;
; PURPOSE:
;
;   This function computes, for each input count energy bin, the fraction of detected counts that result from partial 
;   slat transparency. 
;
; CALLING SEQUENCE:
;
;   corr_factor = stx_slat_transparency_correction(ct_edges, flare_loc, sp_index=sp_index, subc_index=subc_index)
;
; INPUTS:
;
;   ct_edges: float array containing the edges of the count energy bins in keV (e.g., [4, 5, 6, 7] corresponds to 
;             the 4–5 keV, 5–6 keV, and 6–7 keV energy bins).
;   
;   flare_loc: two-element array containing the X and Y coordinates of the estimated flare location (STIX coordinate frame, arcsec). 
;              If not passed, the on-axis grid transmission is computed
;
; OUTPUTS:
;
;   corr_factor: fraction of detected counts attributed to partial slat transparency in each count energy bin.
;   
; KEYWORDS:
;
;   sp_index: float array containing the estimated spectral index of the count spectrum for each energy bin. 
;             By default, the spectral index is equal to 4 in each bin.
;   
;   subc_index:
;
; HISTORY: May 2025, Massa P. and Volpara V., first release
;
; CONTACTS:
;   paolo.massa@fhnw.ch
;   anna.volpara@edu.unige.it
;-

function stx_slat_transparency_correction, ct_edges, flare_loc, sp_index=sp_index, subc_index=subc_index

n_energy_bins = n_elements(ct_edges)-1

default, sp_index, fltarr(n_energy_bins) + 4.
default, subc_index, stx_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c',$
    '6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c'])

;; Define energy edges
energy_low = ct_edges[0:-2]
energy_high = ct_edges[1:-1]

;; Compute slit-to-pitch ratio for front and rear window
subc_transm = stx_subc_transmission(flare_loc)
subc_transm = subc_transm[subc_index]

slit_to_pitch = sqrt(subc_transm)


;; Compute DRM
transmission = read_csv(loc_file( 'stix_transmission_by_component_highres_20240711_010-100eVBin.csv', path = getenv('STX_GRID')))

emin = 1
emax = 150
phe = transmission.field9
phe = phe[where(phe gt emin-1 and phe lt 2*emax)]
edge_products, phe, mean = mean_phe, width = w_phe
;ph_edges = [mean_phe[0] - w_phe[0], mean_phe]
;ph_edges =  get_uniq( [ph_edges,ct_edges],epsilon=0.0001)

;; Read grid thickness
restore,loc_file( 'grid_temp.sav', path = getenv('STX_GRID') )
fff=read_ascii(loc_file( 'grid_param_front.txt', path = getenv('STX_GRID') ),temp=grid_temp)
rrr=read_ascii(loc_file( 'grid_param_rear.txt', path = getenv('STX_GRID') ),temp=grid_temp)

thickness_front_all = fff.THICK
thickness_rear_all = rrr.THICK

thickness_front = fltarr(32)
thickness_rear = fltarr(32)
for subc_n=0,31 do begin

  ;; Front grid
  subc_n_front = fff.SC
  idx = where(subc_n_front eq (subc_n+1), count)
  if count eq 1 then begin

    thickness_front[subc_n] = thickness_front_all[idx]

  endif else begin

    thickness_front[subc_n] = thickness_front_all[idx[0]]

  endelse

  ;; Rear grid
  subc_n_rear = rrr.SC
  idx = where(subc_n_rear eq (subc_n+1), count)
  if count eq 1 then begin

    thickness_rear[subc_n] = thickness_rear_all[idx]

  endif else begin

    thickness_rear[subc_n] = thickness_rear_all[idx[0]]

  endelse

endfor

thickness_front = thickness_front[subc_index]
thickness_rear = thickness_rear[subc_index]

;; Compute energy dependent tranmission
mass_attenuation = xsec(mean_phe, (Element2Z('W'))[0], 'AB', /cm2perg, error=error, use_xcom=1)
gmcm = 19.30
linear_attenuation = mass_attenuation*gmcm/10.

subc_t_energy = fltarr(n_elements(subc_index), n_elements(mean_phe))
for i=0,n_elements(subc_index)-1 do begin
  
  ;; Front grid
  slat_optical_depth = thickness_front[i] * linear_attenuation
  slat_transmission_front = exp(-slat_optical_depth)
  
  ;; Rear grid
  slat_optical_depth = thickness_rear[i] * linear_attenuation
  slat_transmission_rear = exp(-slat_optical_depth)
  
  ;; Subcollimator transmission
  subc_t_energy[i,*]= ( slit_to_pitch[i] * (1-slat_transmission_front) + slat_transmission_front ) * $
                      ( slit_to_pitch[i] * (1-slat_transmission_rear) + slat_transmission_rear ) 
  
endfor


;; Compute approximated correction factor
corr_factor = fltarr(n_elements(subc_index), n_energy_bins)


for i=0,n_energy_bins-1 do begin

  idx = where((mean_phe ge energy_low[i]) and (mean_phe le energy_high[i]))
  this_spect = 1e9 * mean_phe[idx]^(-sp_index[i])

  for j=0,n_elements(subc_index)-1 do begin

      corr_factor[j,i] = f_div(total(reform(subc_t_energy[j,idx]) * this_spect * w_phe[idx]), $
        total(subc_transm[j]* this_spect * w_phe[idx]))

  endfor
endfor

return, corr_factor

end