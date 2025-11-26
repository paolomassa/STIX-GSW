FUNCTION stx_uv_smooth, vis,aux_data,imsize=imsize, pixel=pixel,method=method, threshold_PSI=threshold_PSI, $
                       ep = ep, flare_loc=flare_loc, $
                       NOPLOT = NOPLOT, uv_window = uv_window, _extra =_extra

  ; wrapper around UV_SMOOTH
  ; output map structure has north up
 

  uv_smooth, vis, uv_smooth_im, imsize=imsize, pixel=pixel, method=method, threshold_PSI=threshold_PSI, $
            ep = ep, flare_loc = flare_loc, aux_data = aux_data, $
            NOPLOT = NOPLOT, uv_window = uv_window, _extra =_extra

  this_method        = 'UV_SMOOTH'
  
  
  if (vis[0].type Eq 'stx_visibility') then begin
     uv_smooth_map = stx_make_map(uv_smooth_im.data, aux_data, pixel, this_method, vis)
  endif else begin
     uv_smooth_map = make_map(uv_smooth_im.data)
  endelse
  
  return, uv_smooth_map

END