;+
; NAME:
;   uv_smooth
;
; PURPOSE:
;   This code implements a new interpolation routine for uv-smooth:
;   - utilizes RHESSI or STIX visibilities as input
;   - replaces the unknown visibilities at non-sampled (u,v) points with a VSK interpolant
;   - utilizes an iterative routine incorporating a positivity constraint on the image to provide
;     visibilities that taper gradually to zero outside the sampling domain.
;
; CALLING SEQUENCE:
;   uv_smooth, vis, map
;
; INPUTS:
;   vis:  input visibility structure in standard format
;
; OUTPUTS:
;   map: image map in the structure format provided by the routine make_map.pro
; KEYWORDS:
;   NOPLOT - default to 1 to not plot uv plane of available visibilities
;   imsize: the uv_smooth map size, default to 128X128
;   pixel: the uv_smooth map pixel size, default to [1,1]
;   method: the method used to initialize the procedure (either back-projection or clean component map), 
;     default to 'BP', otherwise 'CC'
;   threshold_PSI: a treshold for the back-projection or clean component map, 
;     default for 'BP' is 0.8, for CC is 0.
;   ep: the shape parameter for the RBF  
;   
; RESTRICTIONS:
;   - input visibilities must be sampled in the whole (u,v) plane (visibilities
;     with only positive v values must be previuosly converted to their conjugate values)
;
; MODIFICATION HISTORY:
;   Dec-2008 Written Anna M. Massone
;   Jul-2009, ras, added RECONSTRUCTED_MAP_VISIBILITIES as an output
;       NOPLOT added to prevent uncontrolled graphics
;       REMOVE_LIST - list of detector isc's not to use because they aren't implemented yet
;   16-Jul-2009, Kim.  Added uv_window arg. Reuse that window for vis sampling plot
;   17-Nov-2009, Anna, Kim. Removed remove_list arg and logic. Can now handle all 9 detectors.
;   26-Nov-2011, Kim. Call al_legend instead of legend (IDL V8 conflict)
;   Dec-2021, Emma Perracchione. Added the VSK interpolation routine and adapted the code to STIX
;-
pro uv_smooth_old, vis, imsize=imsize, pixel=pixel, method=method, threshold_PSI=threshold_PSI, $
  ep = ep, flare_loc = flare_loc, aux_data = aux_data, map, $
  NOPLOT = NOPLOT, uv_window = uv_window, _extra =_extra
;;;;;;;;;;;;;;;;; Set default parameters (almost all inherited by RHESSI) ;;;;;;;;;;;;;;;;;

  default, ep, 0.01
  default, noplot, 1
  default, uv_window, -1
  default, imsize, [128L, 128L]
  default, aux_data, 0
  default, flare_loc, vis[0].xyoffset
  if (vis[0].type Ne 'stx_visibility') and (min(vis.isc) LE 1) then begin    
    default, pixel,  [0.5,0.5]
  endif else begin
    default, pixel, [1.,1.]   
  endelse
  default, method, 'BP'
  if method EQ 'BP' then begin
   default, threshold_PSI, 0.85   
  endif else begin
   default, threshold_PSI, 0.
  endelse
 
;;;;;;;;;;;;;;;;; Set uv-smooth parameters (almost all inherited by RHESSI) ;;;;;;;;;;;;;;;;;
     
  pixel_uv = 0.0005                                 ;;; pixel size in the (u,v)-plane
  imsize_uv_smooth = 128L                           ;;; size of the output image
  Nnew = 1920L                                      ;;; parameter for zero-padded image
  Nmask = 15L                                        ;;; parameter for subsampling the zero-padded image
  iterLand = 50                                     ;;; maximum number of Landweber iterations 
  tau = 0.2                                         ;;; relaxation parameter for Landweber iteration

;;;;;;;;;;;;;;;;; Check if the visibilities are from RHESSI or STIX ;;;;;;;;;;;;;;;;;

  if (vis[0].type Eq 'stx_visibility') then begin
     ;;; Set parameters for STIX visibilities
     fov = 0.16 ;;; set the field of view
     N = 320L   ;;; set the image size (in the (u,v)-plane)
  endif else begin    
    ;;; Set parameters for RHESSI visibilities
    detmin = min(vis.isc)  
    if (detmin eq 0 ) then begin
      fov = 0.45 ;;; set the field of view
      pixel_uv = pixel_uv*2.
      N = 450L   ;;; set the image size (in the (u,v)-plane)
    endif
    if (detmin eq 1 ) then begin
      fov = 0.26 ;;; set the field of view
      pixel_uv = pixel_uv*2.
      N = 260L   ;;; set the image size (in the (u,v)-plane)
    endif
    if (detmin GE 2 ) then begin
      fov = 0.16 ;;; set the field of view
      N = 320L   ;;; set the image size (in the (u,v)-plane)
      endif
  endelse

;;;;;;;;;;;;;;;;; Plot of the visibility sampling on the (u,v)-plane ;;;;;;;;;;;;;;;;;
  
  if not NOPLOT then begin
    ;;; Save current window. Reuse old uv_window if available, otherwise create it.
    save_window=!d.window
    if uv_window eq -1 or (is_wopen(uv_window) eq 0) then begin
      uv_window = next_window(/user)
      window, uv_window, xsize=500, ysize=500,xpos=0,ypos=50,title='u-v sampling '
    endif
    wset, uv_window
    plot, vis.u, vis.v,  /isotropic, psym=8, xtitle='u (arcsec!u-1!n)', ytitle='v (arcsec!u-1!n)', $
      title='Visibility Sampling in u-v plane'
    if (vis[0].type Eq 'stx_visibility') then begin
      al_legend, [format_intervals(minmax(anytim(vis[0].time_range[0].value.time, $
        mjd = vis[0].time_range[0].value.mjd, /ecs)),/ut)], box=0
    endif else begin  
      det = 'Detectors: ' + arr2str(trim(get_uniq(vis.isc)+1), ' ')
      al_legend, [format_intervals(minmax(vis.trange),/ut), det], box=0
    endelse
    ;;; Reset to user's window
    wset,save_window
  endif
  
  ;;;;;;;;;;;;;;;;; Definition of the uniform grid for the interpolation ;;;;;;;;;;;;;;;;;

  Ulimit = (N/2.-1)*pixel_uv+pixel_uv/2.
  usampl = -Ulimit+findgen(N)*pixel_uv
  vsampl = usampl

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;; Interpolation with VSKs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

reintz= GRID_TPS (vis.u, vis.v, float(vis.obsvis)/(4*!pi*!pi), NGRID=[N,N] , START=[usampl[0], vsampl[0]], DELTA=[pixel_uv,pixel_uv] )
imintz=GRID_TPS (vis.u, vis.v, imaginary(vis.obsvis)/(4*!pi*!pi), NGRID=[N,N] , START=[usampl[0], vsampl[0]], DELTA=[pixel_uv,pixel_uv] )
visnew=complex(reintz,imintz)


  ;;;;;;;;;;;;;;;;; Do not allow function to extrapolate outside the disk  ;;;;;;;;;;;;;;;;;
  
  Rmax = max(sqrt(vis.u*vis.u+vis.v*vis.v))
  for i = 0L,N-1 do begin
    for j = 0L,N-1 do begin
      if(sqrt(usampl[i]*usampl[i]+vsampl[j]*vsampl[j]) GT Rmax ) then visnew[i,j] = complex(0.,0.)
    end
  end

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;; Zero-padding and subsampling
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;; First step: Zero-padding  ;;;;;;;;;;;;;;;;;

  Ulimit = (Nnew/2.-1)*pixel_uv+pixel_uv/2.
  xpix = -Ulimit+findgen(Nnew)*pixel_uv
  ypix = xpix
  intzpadd = MAKE_ARRAY(Nnew,Nnew, /complex, VALUE = 0)
  intzpadd[(Nnew-N)/2.:(Nnew-N)/2.+N-1,(Nnew-N)/2.:(Nnew-N)/2.+N-1] = visnew
    
  ;;;;;;;;;;;;;;;;;  Second step: resampling of the image with a mask ;;;;;;;;;;;;;;;;
  
  im_new = Nnew/Nmask
  intznew = complexarr(im_new,im_new)
  xpixnew = fltarr(im_new)
  ypixnew = fltarr(im_new)

  for i = 0,im_new-1 do begin
    xpixnew[i] = xpix[Nmask*i+Nmask/2]
    ypixnew[i] = ypix[Nmask*i+Nmask/2]
    for j=0,im_new-1 do intznew[i,j] =  intzpadd[Nmask*i+Nmask/2,Nmask*j+Nmask/2]
  end

  ;;;;;;;;;;;;;;;;; Arrange parameter values before applying the FFT  ;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;; Compute the fov and the sampling distance ;;;;;;;;;;;;;;;;;
  
  OMEGA = (xpixnew[im_new-1]-xpixnew[0])/2.
  X = im_new/(4*OMEGA)
  deltaomega = (xpixnew[im_new-1]-xpixnew[0])/im_new
  deltax = 2*X/im_new

  ;;;;;;;;;;;;;;;;; Characteristic function of the disk ;;;;;;;;;;;;;;;;;

  chi = complexarr(im_new,im_new)
  for i = 0,im_new-1 do begin
    for j = 0,im_new-1 do begin
      if (sqrt(xpixnew[i]*xpixnew[i]+ypixnew[j]*ypixnew[j]) LE Rmax) then chi[i,j] = complex(1.,0.)
    end
  end
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;; Projected Landweber method with positivity constraint
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  map_actual = complexarr(im_new,im_new)
  map_iteration = fltarr(iterLand,im_new,im_new)
  map_solution = fltarr(im_new,im_new)
  descent = fltarr(iterLand-1)
  normAf_g = fltarr(iterLand)

  ;;;;;;;;;;;;;;; iteration 0: Inverse Fourier Transform of the initial solution ;;;;;;;;;;;;;;;;;
  
  map_shifted = shift(map_actual,[-im_new/2.,-im_new/2.])
  F_Trasf_shifted = FFT(map_shifted,/inverse)/(4*!pi*!pi*deltaomega*deltaomega*im_new*im_new)
  F_Trasf = shift(F_Trasf_shifted,[im_new/2.,im_new/2.])

  ;;;;;;;;;;;;;;; Landweber i-th iteration ;;;;;;;;;;;;;;;;;
  
  for iter = 0,iterLand-1 do begin
    F_Trasf_up = F_Trasf+ tau*(intznew- chi*F_Trasf)          ;;;; Landweber updating rule
    F_Trasf = F_Trasf_up
    F_Trasf_shifted = shift(F_Trasf,[-im_new/2.,-im_new/2.])  ;;;; Fourier Transform of the updated solution
    map_shifted = FFT(F_Trasf_shifted)*4*!pi*!pi*deltaomega*deltaomega*im_new*im_new
    map_actual = shift(map_shifted,[im_new/2.,im_new/2.])
    ;;;; Projection of the solution onto the subset of the positive solutions (positivity constraint)
    for i = 0,im_new-1 do begin
      minzero = where( float(map_actual[0:im_new-1,i]) lt 0,count)
      if (count NE 0) then  map_actual[minzero,i] = complex(0.,0.)
    endfor
    map_iteration[iter,*,*] = float(map_actual)
    map_shifted = shift(map_actual,[-im_new/2.,-im_new/2.])
    F_Trasf_shifted = FFT(map_shifted,/inverse)/(4*!pi*!pi*deltaomega*deltaomega*im_new*im_new)
    F_Trasf = shift(F_Trasf_shifted,[im_new/2.,im_new/2.])
    ;;;;;;;;;;;;;;;; Stop criterion based on the descent of ||Af-intznew||
    Af_g = chi*F_Trasf-intznew
    normAf_g[iter] = sqrt(total(abs(Af_g)*abs(Af_g)))
    if (iter GE 1) then begin
      descent[iter-1] = (normAf_g[iter-1]-normAf_g[iter])/normAf_g[iter-1]
      if (descent[iter-1] LT 0.02  ) then break
    endif
  endfor

  ;;;;;;;;;;;;;;;;;;; Map corresponding to the optimal iteration ;;;;;;;;;;;;;;;;;
  
  if (iter eq iterLand) then map_solution[*,*] = map_iteration[14,*,*]  else map_solution = float(map_actual)

  ;;;;;;;;;;;;;;;;; Consistency check for RHESSI and STIX time visibilities ;;;;;;;;;;;;;;;;; 
  
  flag = 1
  catch, error_status
  if error_status ne 0 then begin
    aux = vis[0].time_range[0].value
    time = anytim(aux.time, mjd = aux.mjd, /ecs)
    flag = 0
    catch, /cancel
  endif
  if flag then time = anytim(vis[0].trange[0], /ecs)

  ;;;;;;;;;;;;;;;;;;; Make the map ;;;;;;;;;;;;;;;;;

  mapcenter = flare_loc
  map = make_map(map_solution,$
    id=' ', $ ;earth orbit
    xc=mapcenter[0],yc=mapcenter[1],$
    dx=deltax,dy=deltax,$
    time=time, $
    xunits='arcsec', yunits='arcsec')
  
  ;;;;;;;;;;;;;;;;;;; If needed rebin the map ;;;;;;;;;;;;;;;;;
  
  if (vis[0].type Ne 'stx_visibility') and (min(vis.isc) LE 1) then begin    
    if same_data(fix(imsize), [128,128]) NE 1 $
      or same_data(pixel, [0.5,0.5]) NE 1 then begin
       newmap = rep_tag_value(map, fltarr(imsize),'DATA')
      newmap.dx = pixel[0] & newmap.dy = pixel[1]
      map = inter_map(map, newmap, err=err)
    endif  
  endif else begin
  if same_data(fix(imsize), [128,128]) NE 1 $
     or same_data(pixel, [1.,1.]) NE 1 then begin
    newmap = rep_tag_value(map, fltarr(imsize),'DATA')
    newmap.dx = pixel[0] & newmap.dy = pixel[1]
    map = inter_map(map, newmap, err=err)
  endif
  endelse
end


