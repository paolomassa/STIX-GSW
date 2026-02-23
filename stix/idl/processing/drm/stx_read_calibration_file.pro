;+
;
; name:
;       stx_read_calibration_file
;
; :description:
;    Read the values contained in a STIX calibration file and returns a correponding 'stx_calibration_data' structure
;    
;
; :params:
;    fits_path : in, required, type="string"
;                the path of the FITS file to be read. Passed through to mrdfits.
;
;
; :keywords:
;
;    primary_header : out, type="string array"
;              The header of the primary HDU of the pixel data files
;
;    data_header : out, type="string array"
;              The header of the data extension of the auxiliary file
;
;    data_str : out, type="structure"
;              The contents of the data extension of the auxiliary file
;
;    control_header : out, type="string array"
;                The header of the control extension of the auxiliary file
;
;    control_str : out, type="structure"
;               The contents of the control extension of the auxiliary file
;               
;    idb_version_header : out, type="string array"
;               The header of the idb version extension of the auxiliary file
;               
;    idb_version_str : out, type="structure"
;               The contents of the idb version extension of the auxiliary file
;
; :returns:
;
;    a 'stx_calibration_data' structure containing the values read from a STIX calibration fits file
;
; :examples:
;
;    calib_data = stx_read_aux_fits(fits_path)
;
; :history:
;
;    August 2025: Massa P. (FHNW), created
;
;-

function stx_read_calibration_file, fits_path, primary_header = primary_header, data_str = data, $
                                    data_header = data_header, control_str = control, control_header= control_header, $
                                    idb_version_str = idb_version, idb_version_header = idb_version_header
  
  
  ;; Read fits file
  !null       = stx_read_fits(fits_path, 0, primary_header,  mversion_full = mversion_full, /silent)
  control     = stx_read_fits(fits_path, 'control', control_header, mversion_full = mversion_full, /silent)
  data        = stx_read_fits(fits_path, 'data', data_header, mversion_full = mversion_full, /silent)
  idb_version = stx_read_fits(fits_path, 'idb_versions', idb_version_header, mversion_full = mversion_full, /silent)

  ;; Create time object
  start_time = sxpar(primary_header, 'date-beg')
  stx_time_obj = stx_time()
  stx_time_obj.value =  anytim(start_time , /mjd)
  
  time_bin_center = data.TIME ;; in cs
  duration = data.TIMEDEL ;; in cs
  
  t_start = stx_time_add( stx_time_obj,  seconds = [ time_bin_center/100 - duration/200 ] )
  t_end   = stx_time_add( stx_time_obj,  seconds = [ time_bin_center/100 + duration/200 ] )
  t_mean  = stx_time_add( stx_time_obj,  seconds = [ time_bin_center/100 ] )

  ;; Create calibration structure
  calib_data = stx_calibration_data()
  
  calib_data.time_start = t_start
  calib_data.time_end = t_end
  calib_data.t_mean = t_mean
  
  calib_data.energy_bin_low = data.E_EDGES_ACTUAL[0:31,*,*]
  calib_data.energy_bin_high = data.E_EDGES_ACTUAL[1:32,*,*]
  
  calib_data.gain = 1./reform(data.GAIN)
  calib_data.offset = reform(data.OFFSET)
  
  calib_data.live_time = data.LIVE_TIME
  
  calib_data.elut_name = control.OB_ELUT_NAME
  
  return, calib_data

end
