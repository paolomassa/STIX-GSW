;+
;
; NAME:
;   distance_matrix
;
; PURPOSE:
;   Computes the distance matrix between two sets of points
;
; CALLING SEQUENCE:
;   distance_matrix, dsites, ctrs
;
; INPUTS:
;   dsites:    a set of points
;   ctrs:      a set of RBF centres
;
; OUTPUTS:
;   DM:  the distance matrix
;
; NOTES;
;   The function is taken by [G.E. Fasshauer, Meshfree Approximation Methods 
;   with Matlab, World Scientific, Singapore, 2007]
;
; HISTORY: Dec 2021, Perracchione E., created
;
function distance_matrix, dsites, ctrs
  
;;;;;;;;;;;;;;;;; Initialize ;;;;;;;;;;;;;;;;;
  NN =  SIZE(dsites[*,0], /N_ELEMENTS)
  NN1 =  SIZE(ctrs[*,0], /N_ELEMENTS)
  s = SIZE(dsites[0,*], /N_ELEMENTS)
  DM = FLTARR(NN1, NN)
  
  ;;;;;;;;;;;;;;;;; Compute the distance matrix ;;;;;;;;;;;;;;;;;
  for d=0L,s-1 do begin
    dr =  FLTARR(NN, NN1)
    cc =  FLTARR(NN, NN1)
    dr = transpose(replicate_array(dsites[*,d], NN1))
    cc = replicate_array(ctrs[*,d], NN)
    ;  dr, consisting of N identical columns (each containing
    ;      the d-th coordinate of the M data sites)
    ;  cc, consisting of M identical rows (each containing
    ;      the d-th coordinate of the N centers)
    DM = DM + (dr-cc)^2
  end
  
  DM = sqrt(DM) 
  return, DM

end