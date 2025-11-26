;+
;
; NAME:
;   matern_kernel_interp
;
; PURPOSE:
;   Kernel interpolation with the Matern RBF
;
; CALLING SEQUENCE:
;   matern_kernel_interp, dsites, epoints, rhs, ep, regparam, multi_set
;
; INPUTS:
;   dsites:    a set of points
;   epoints:   a set of RBF centres
;   rhs:       a set of data values (e.g. visibilities)
;   ep:        the shape parameter of the RBF interpolant
;   regparam:  regression parameter
;   multi_set: useful if the interpolant has to be evaluated on multiple datasets
;
; OUTPUTS:
;   Pf:  the interpolant
;
; NOTES;
;   Part of the function is taken by [G.E. Fasshauer, Meshfree Approximation Methods
;   with Matlab, World Scientific, Singapore, 2007]
;
; HISTORY: Dec 2021, Perracchione E., created
;
function matern_kernel_interp, dsites, epoints, rhs, ep, regparam = regparam, $
  multi_set = multi_set

;;;;;;;;;;;;;;;;; Set the parameters ;;;;;;;;;;;;;;;;;

  default, multi_set,  0
  default, regparam, fltarr(size(dsites[*,0],/N_ELEMENTS)) 

  ;;;;;;;;;;;;;;;;; Compute the kernel matrix for interpolation ;;;;;;;;;;;;;;;;;

  DM_kernel = distance_matrix(dsites, dsites)  
  IM = exp(-ep*DM_kernel)
   
  ;;;;;;;;;;;;;;;;; Solve the interpolation system ;;;;;;;;;;;;;;;;;
  
  coef = fltarr(size(rhs,/dim))
  for i=0L,size(rhs[0,*], /N_elements)-1 do begin
  coef[*,i] = reform(LA_LINEAR_EQUATION(IM+diag_matrix(transpose(regparam)),$
    transpose(rhs[*,i])),size(dsites[*,0],/N_ELEMENTS),1)
  endfor

  ;;;;;;;;;;;;;;;;; Compute the evaluation matrix for interpolation ;;;;;;;;;;;;;;;;;
  
  NN1 =  SIZE(epoints[*,0], /N_ELEMENTS)
  N_max = 450^2./2
  T = ceil(NN1/N_max)
  Pf = fltarr(NN1,size(rhs[0,*], /N_elements))
  count = 0.
  k = 1.
  while k LE T do begin 
    idx = min([NN1,N_max])
    if k*idx LE NN1 then begin
      idx_1 = k*idx
    endif else begin
      idx_1 = NN1
    endelse
    DM = distance_matrix(dsites, epoints[count:idx_1-1,*])  
    EM = exp(-ep*DM)
    for nsys=0,size(rhs[0,*], /N_elements)-1 do begin
      Pf[count:idx_1-1,nsys] = (EM)#coef[*,nsys]
    endfor
          
  count += idx
  k += 1
endwhile
          
if max(multi_set) NE 0 then begin  
  NN1 =  SIZE(multi_set[*,0], /N_ELEMENTS)
  T = ceil(NN1/N_max)
  Pf_2 = fltarr(NN1,size(rhs[0,*], /N_elements))
  count = 0.
  k = 1.
  while k LE T do begin 
    idx = min([NN1,N_max])
    if k*idx LE NN1 then begin
      idx_1 = k*idx
    endif else begin
      idx_1 = NN1
    endelse
    DM = distance_matrix(dsites, multi_set[count:idx_1-1,*])  
    EM = exp(-ep*DM)
    for i=0L,size(rhs[0,*], /N_elements)-1 do begin
      Pf_2[count:idx_1-1,i] = (EM)#coef[*,i]
    endfor
  count += idx
  k += 1
endwhile

      
 if max(multi_set) NE 0 then begin
    Pf = {Pf1:Pf, Pf2:Pf_2}
  endif
endif

return, Pf   

   end