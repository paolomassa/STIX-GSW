;+
;
; NAME:
;   uv_smooth_VSK
;
; PURPOSE:
;   This code computes a VSK interpolant based on the usage of a first rought image reconstruction
;
; CALLING SEQUENCE:
;   uv_smooth_vsk, vis, usampl, N, method, Ulimit, perc_PSI, imsize_default, Nnew, fov, pixel_uv
;
; INPUTS:
;   vis: the visibilities
;   usampl: the evaluation data (i.e. the vector defining the grid)
;   N: number of grid data in one direction
;   method: the method used to initialize the procedure (either back-projection or clean component map)
;   Ulimit: the range of the square cointaining the data in the uv-plane
;   perc_PSI: a treshold for the back-projection or clean component map
;   imsize_default: the uv_smooth map size
;   Nnew: fov after the zero padding strategy
;   fov: fov in the uv-plane
;
; OUTPUTS:
;   visnew: the visibility surfaces
;
; NOTES: Part of this code is taken by some Matlab codes available at https://github.com/emmaA89/VSDKs/
;
function normalize_data, data, data1 = data1, second_set = second_set

  if KEYWORD_SET(data1) EQ 0 then begin
     data1 = data
     if KEYWORD_SET(second_set) then begin
       data_rescaled = data*max(second_set-min(second_set)) + min(second_set)
     endif else begin
       data_rescaled = (data-min(data1))/max(data1-min(data1))
     endelse
  endif else begin
  if KEYWORD_SET(second_set) then begin
    data_rescaled = data*max(second_set-min(second_set)) + min(second_set)
  endif else begin
    data_rescaled = (data-min(data1))/max(data1-min(data1))
  endelse
  endelse

  return, data_rescaled

end
