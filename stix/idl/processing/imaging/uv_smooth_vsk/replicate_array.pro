;+
;
; NAME:
;   replicate_array
;
; PURPOSE:
;   Replicate an array n times
;
; CALLING SEQUENCE:
;   replicate_array, arr, n
;
; INPUTS:
;   arr:    an array
;   n:      number of times to replicate
;
; OUTPUTS:
;   replicated_array:  array consisting of arr replicated n times
;
;

function replicate_array, arr, n

dim = size(arr, /dimensions)
if dim eq 0 then begin
  arr1 = replicate(arr, n)  
endif else begin

replicated_array = fltarr(dim, n)

for i=0, n-1 do begin

replicated_array[*, i] = arr  
  
endfor
endelse
return, replicated_array

end