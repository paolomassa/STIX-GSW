;+
;
; name:
;       stx_get_calibration_file
;
; :description:
;    This procedure checks the STIX data archive for a given observation time and finds any
;    STIX CAL calibration file. If no file is present for the input date, the procedure searches for 
;    the calibration file which is closest in time. If multiple files are present for that date,
;    this procedure downloads the one with largest duration. 
;
; :categories:
;    template, example
;
; :params:
;    start_time : in, required, type="string"
;                the start time of the observation
;    end_time : in, required, type="string"
;                the end time of the observation
;
; :keywords
;   out_dir: path of the folder where the STIX calibration FITS files are saved. Default is the current directory
;
;   clobber: 0 or 1. If set to 0, the code does not download the file again if it is already present in 'out_dir'.
;
; :returns:
;   Path of the downloaded STIX calibration FITS file
;
; :examples:
;    out_file = stx_get_calibration_file('09-May-23 06:14:37.094', '09-May-23 06:36:12.194')
;
; :history:
;    26-Jan-2026 - Massa P. (FHNW), first release
;-
function stx_get_calibration_file, start_time, end_time, out_dir=out_dir, clobber=clobber

  cd, current=current

  default, out_dir, current
  default, clobber, 0

  site = 'http://dataarchive.stix.i4ds.net'
  date_path = get_fid(start_time,end_time,/full,delim='/')

  type_path = '/testing/p311_comp_test/CAL/'
  path  = type_path + date_path[0] +'/CAL'
  filter = '*stix-cal-energy*.fits'
  found_files=sock_find(site,filter,path=path,count=count)
  
  if count eq 0 then begin
    
    message, "No STIX calibration files found for this day ("+date_path[0]+"). Looking for the closest in time.", /continue
    
    duration_day = 86400.0 ;; in sec
   
    ;; Back in time
    i_back = 1
    count_back = 0
    while (count_back eq 0) and (i_back le 30) do begin  ;; Go back for a maximum of 30 days
      
      this_start_time = anytim(anytim(start_time) - i_back * duration_day,/vms)
      this_end_time = anytim(anytim(end_time) - i_back * duration_day,/vms)
      
      date_path_back = get_fid(this_start_time,this_end_time,/full,delim='/')
      path  = type_path + date_path_back[0] +'/CAL'
      found_files_back=sock_find(site,filter,path=path,count=count_back)
      
      i_back += 1
      
    endwhile
    
    ;; Forward in time
    i_forward = 1
    count_forward = 0
    while (count_forward eq 0) and (i_forward le 30) do begin ;; Go forward for a maximum of 30 days
      
      this_start_time = anytim(anytim(start_time) + i_forward * duration_day,/vms)
      this_end_time = anytim(anytim(end_time) + i_forward * duration_day,/vms)
      
      date_path_forward = get_fid(this_start_time,this_end_time,/full,delim='/')
      path  = type_path + date_path_forward[0] +'/CAL'
      found_files_forward=sock_find(site,filter,path=path,count=count_forward)
      
      i_forward += 1
      
    endwhile
    
    ;; No file found
    if (count_back eq 0) and (count_forward eq 0) then message, "No STIX calibration file found for this day ("+date_path[0]+")"
      
    if (count_back gt 0) and (count_forward eq 0) then begin
      
      found_files = found_files_back
      count = count_back
      date_path = date_path_back
      
    endif
    
    if (count_back eq 0) and (count_forward gt 0) then begin

      found_files = found_files_back
      count = count_back
      date_path = date_path_back

    endif
      
    if (count_back gt 0) and (count_forward gt 0) then begin
      
      if count_back lt count_forward then begin
        
        found_files = found_files_back
        count = count_back
        date_path = date_path_back
        
      endif else begin
        
        found_files = found_files_back
        count = count_back
        date_path = date_path_back
        
      endelse
  
    endif
  
  endif
  
  message, "Download STIX calibration file recorded at " +date_path[0], /continue
    
  if count eq 1 then begin
  
  selected_file = found_files[0]
  sock_copy, selected_file, out_name, local_file=out_file, out_dir = out_dir, clobber=clobber
  
  endif else begin
    
    message, "More than 1 STIX calibration file found for this day ("+date_path[0]+"). Looking for the one with largest duration.", /continue
    
    start_time_file = []
    end_time_file = []
    
    ;; Extract file names
    for i = 0,n_elements(found_files)-1 do begin
      len_path = STRLEN(site+path)
      len_full_path = STRLEN(found_files[i])
      filename = STRMID(found_files[i], len_path+1, len_full_path)
      
      string_dates = STRMID(filename, STRLEN('solo_CAL_stix-cal-energy_'), STRLEN(filename)-STRLEN('solo_CAL_stix-cal-energy_')-STRLEN('_V02.fits'))
      
      mid_part = STRPOS(string_dates, '-')
      
      this_start_time_file = STRMID(string_dates, 0, mid_part)
      year  = STRMID(this_start_time_file, 0, 4)
      month = STRMID(this_start_time_file, 4, 2)
      day   = STRMID(this_start_time_file, 6, 2)
      hour  = STRMID(this_start_time_file, 9, 2)
      min   = STRMID(this_start_time_file, 11, 2)
      sec   = STRMID(this_start_time_file, 13, 2)
      start_time_file = [start_time_file, anytim(year + '-' + month + '-' + day + 'T' + hour + ':' + min + ':' + sec)]
      
      this_end_time_file = STRMID(string_dates, mid_part+1, STRLEN(string_dates))
      year  = STRMID(this_end_time_file, 0, 4)
      month = STRMID(this_end_time_file, 4, 2)
      day   = STRMID(this_end_time_file, 6, 2)
      hour  = STRMID(this_end_time_file, 9, 2)
      min   = STRMID(this_end_time_file, 11, 2)
      sec   = STRMID(this_end_time_file, 13, 2)
      end_time_file = [end_time_file, anytim(year + '-' + month + '-' + day + 'T' + hour + ':' + min + ':' + sec)]
    
    endfor
    
    ;; Select file with largest duration
    duration = end_time_file - start_time_file
    
    idx_file = where(duration eq max(duration))
    
    selected_file = found_files[idx_file]
    sock_copy, selected_file, out_name, local_file=out_file, out_dir = out_dir, clobber=clobber
    
    
  endelse
    
  return, out_file
  
end