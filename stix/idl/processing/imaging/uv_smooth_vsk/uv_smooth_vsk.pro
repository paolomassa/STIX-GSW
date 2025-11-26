;+
;
; NAME:
;   uv_smooth_VSK
;
; PURPOSE:
;   This code computes a VSK interpolant based on the usage of a first rought image reconstruction
;
; CALLING SEQUENCE:
;   uv_smooth_vsk, vis, usampl, N, method, Ulimit, threshold_PSI, ep, aux_data 
;   
; INPUTS:
;   vis: the visibilities 
;   usampl: the evaluation data (i.e. the vector defining the grid)
;   N: number of grid data in one direction
;   method: the method used to initialize the procedure (either back-projection or clean component map)
;   Ulimit: the range of the square cointaining the data in the uv-plane
;   threshold_PSI: a treshold for the back-projection or clean component map
;   ep: the shape parameter for the RBF  
;   aux_data: aux file of the visibilities (mandatory for clean map)
;
; OUTPUTS:
;   visnew: the visibility surfaces
; 
; NOTES: Part of this code is taken by some Matlab codes available at https://github.com/emmaA89/VSDKs/
;
; HISTORY: Dec 2021, Perracchione E., created
;

function uv_smooth_vsk, vis, usampl, N, method, Ulimit, threshold_PSI, ep, aux_data 

;;;;;;;;;;;;;;;;; Store the data in matrices and vectors for interpolation ;;;;;;;;;;;;;;;;;

  vis_matrix = [[(vis.u)], [(vis.v)]]   
  vis_real = normalize_data(float(vis.obsvis))
  vis_imaginary = normalize_data(imaginary(vis.obsvis))
  
;;;;;;;;;;;;;;;;; Define the evaluation points: a matrix (N^2 X 2) ;;;;;;;;;;;;;;;;;

  grid_matrix = [[reform(transpose(replicate_array(usampl, N)),N^2,1)], $
    [reform(replicate_array(usampl, N),N^2,1)]]
     
;;;;;;;;;;;;;;;;; Define the regression parameter ;;;;;;;;;;;;;;;;; 
;;;;;;;;;;;;;;;;; It is non-zero only for RHESSI detectors 1 and 2 ;;;;;;;;;;;;;;;;; 
  
  Nall = size(vis_matrix[*,0],/N_ELEMENTS)
  regparam = reform(transpose(0*findgen(Nall[0])),Nall[0])
  wh = [0:size(vis_matrix[*,0],/N_ELEMENTS)-1]
  if (vis[0].type Ne 'stx_visibility') then begin
      if min(vis.isc) LE 1 then begin
         snr_vis = hsi_vis_get_snr(vis)
         wh = where(vis.isc GE 2)
         lambda_reg = 0.01/snr_vis
         regparam[where(vis.isc LE 1)] = lambda_reg
      endif
  endif

;;;;;;;;;;;;;;;;; Compute a first approximation of the image ;;;;;;;;;;;;;;;;; 

   if (method eq 'BP') then begin
           if (vis[0].type Ne 'stx_visibility') then begin
               vis_bpmap, vis[wh], MAP=map, data_only = 0
               dx = map.dy
               dy = map.dx
               map_coarse =   map.data
           endif else begin
               bp_map = stx_bproj(vis[wh], [128., 128.], [1., 1.], aux_data)
               a = rot_map(bp_map,aux_data.ROLL_ANGLE,rcenter=[0.,0.])
               map_coarse =   rotate(a.data,1)
               px_1 = bp_map.dx
               dx = px_1[0]
               dy =  dx
             endelse
   endif else begin 
           if (vis[0].type Ne 'stx_visibility') then begin
              clean_im = vis_clean(vis[wh], clean_map = clean_map, $
              clean_sources_map = clean_sources_map, pixel=px_1)
              map_coarse =   clean_sources_map
              dx = px_1[0]
              dy =  dx
           endif else begin
             clean_im = stx_vis_clean(vis[wh], aux_data, clean_map = clean_map, $
             clean_sources_map = clean_sources_map)
             a = rot_map(make_map(clean_sources_map),aux_data.ROLL_ANGLE,rcenter=[0.,0.])
             map_coarse =   rotate(a.data,1)
             px_1 = clean_im.dx
             dx = px_1[0]
             dy =  dx
           endelse
    endelse

;;;;;;;;;;;;;;;;; Threshold the rought map ;;;;;;;;;;;;;;;;; 
   N_m = sqrt(size(map_coarse,/N_elements))
   tres = MAKE_ARRAY(N_m, N_m, /COMPLEX, VALUE = threshold_PSI*max(map_coarse))
   map_coarse_tres =  (map_coarse)*(abs(map_coarse) ge tres)
   if size(where(map_coarse_tres GT 0),/N_elements) EQ 1 then map_coarse_tres = map_coarse
   
;;;;;;;;;;;;;;;;; Compute the augmented features ;;;;;;;;;;;;;;;;; 
   augemted_features = uv_smooth_augmented_feature(map_coarse_tres, Ulimit, grid_matrix, $
    vis_matrix, ep, N_m, dx, dy, method, aux_data , vis[wh])
  
;;;;;;;;;;;;;;;;; Interpolate the visibilities ;;;;;;;;;;;;;;;;;    
   interp_vis = matern_kernel_interp(augemted_features.vis_matrix_VSK, augemted_features.grid_matrix_VSK, [[vis_real], [vis_imaginary]], ep, $
     regparam = regparam, multi_set = multi_set)

;;;;;;;;;;;;;;;;; Compute the complex grids ;;;;;;;;;;;;;;;;  
   reintz = reform(normalize_data(interp_vis[*,0], second_set = float(vis.obsvis))/$
     (4*!pi*!pi), N, N)  
   imintz = reform(normalize_data(interp_vis[*,1], second_set = imaginary(vis.obsvis))/$
     (4*!pi*!pi), N, N)  
   visnew = rotate(complex(reintz,imintz),4)
    
   return, visnew
  
  end
  
 