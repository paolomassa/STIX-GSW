;+
;
; NAME:
;   uv_smooth_augmented_feature
;
; PURPOSE:
;   This code computes the added features 
;
; CALLING SEQUENCE:
;   uv_smooth_augmented_feature, map_coarse_tres, Ulimit, grid_matrix, vis_matrix, ep,  $
;   N_m, dx, dy, method, aux_data, vis_wh
;   
; INPUTS:
;   map_coarse_tres: the thresholded rought map
;   Ulimit: the range of the square cointaining the data in the uv-plane
;   grid_matrix: the grid data 
;   vis_matrix: the visibilities
;   ep: the shape parameter for the RBF     
;   N_m: size of map_coarse_tres
;   dx, dy: pixel size of the map_coarse_tres
;   method: the method used to initialize the procedure (either back-projection or clean component map),
;     default to 'BP', otherwise 'CC'
;   aux_data: aux file of the visibilities (mandatory for clean map)
;   vis_wh: usefull for RHESSI and is related to smoothing parameter
;   
; OUTPUTS:
;   extra_features: the added features for both visibilities and grid points
;   
; HISTORY: Dec 2021, Perracchione E., created
;

function uv_smooth_augmented_feature, map_coarse_tres, Ulimit, grid_matrix, vis_matrix, ep, $
   N_m, dx, dy, method, aux_data, vis_wh
   
;;; Use the map vis_map2vis_matrix to have a first approximation of the visibility surfaces   
   N_tmp = 24 
   grid_matrix_N_tmp = [[reform(transpose(replicate_array(-Ulimit+findgen(N_tmp)*(2*Ulimit/N_tmp), N_tmp)),N_tmp^2,1)],$
     [reform(replicate_array(-Ulimit+findgen(N_tmp)*(2*Ulimit/N_tmp), N_tmp),N_tmp^2,1)]]
   ; Compute the augmented feature for the data sites
   F = vis_map2vis_matrix([grid_matrix_N_tmp[*,0]], [grid_matrix_N_tmp[*,1]], $
    [N_m, N_m], [dx, dy])
   grid_matrix_tmp_augmented = F # map_coarse_tres[*]

;;;;;;;;;;;;;;;;; Normalize data before interpolation ;;;;;;;;;;;;;;;;;
   grid_matrix_tmp_augmented_real = normalize_data(float(grid_matrix_tmp_augmented)) 
   grid_matrix_tmp_augmented_imag = normalize_data(imaginary(grid_matrix_tmp_augmented)) 
    
 ;;;;;;;;;;;;;;;;; Compute the augmented features ;;;;;;;;;;;;;;;;;   
   augmented_features = matern_kernel_interp(grid_matrix_N_tmp, grid_matrix, $
    [[grid_matrix_tmp_augmented_real], [grid_matrix_tmp_augmented_imag]], ep, $
     regparam = regparam, multi_set = vis_matrix)

 ;;;;;;;;;;;;;;;;; Store the features ;;;;;;;;;;;;;;;;;  
   grid_matrix_VSK = [[grid_matrix] ,[normalize_data([augmented_features.Pf1], $
    data1=[augmented_features.Pf1])]]
   vis_matrix_VSK = [[vis_matrix],  [normalize_data([augmented_features.Pf2], $
    data1=[augmented_features.Pf1])]]

   extra_features = {name:'', vis_matrix_vsk:vis_matrix_VSK, grid_matrix_vsk:grid_matrix_VSK} 

return, extra_features

end