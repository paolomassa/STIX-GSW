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
;    24-Mar-2026 - Massa P. (FHNW), first release
;    13-Apr-2026 - Massa P. (FHNW), updated to search for calibration files within 2 days before or after the input date
;-
function stx_get_calibration_file, start_time, end_time, out_dir=out_dir, clobber=clobber

  cd, current=current

  default, out_dir, current
  default, clobber, 0

  day_in_s = 86400.d ;; Number of seconds in a day. To be used later to identify the most appropriate calibration file
  day_limit = 2 ;; Maximum number of days before or after the input date to consider when searching for a calibration file (arbitrary choice)

  site = 'http://dataarchive.stix.i4ds.net'
  date_path = get_fid(start_time,end_time,/full,delim='/')

  type_path = '/fits/CAL/'
  path  = type_path + date_path[0] +'/CAL'
  filter = '*stix-cal-energy*.fits'
  found_files=sock_find(site,filter,path=path,count=count)

  if count eq 0 then begin
    message, [" ", " ", "No STIX calibration files found for this day ("+date_path[0]+"). Searching for calibration file which is closest in time.", " ", " "], /continue
    
    ;; Forward in time (for 2 days max - arbitrary choice)
    day_n_forward=1
    count_forward=0
    while (day_n_forward le day_limit) and (count_forward eq 0) do begin
      
      this_start_time_forward = anytim(anytim(start_time) + day_n_forward * day_in_s, /vms)
      this_end_time_forward   = anytim(anytim(end_time) + day_n_forward * day_in_s, /vms)
      
      date_path_forward = get_fid(this_start_time_forward,this_end_time_forward,/full,delim='/')
      
      path  = type_path + date_path_forward[0] +'/CAL'
      found_files_forward = sock_find(site,filter,path=path,count=count_forward)
      
      day_n_forward += 1
      
    endwhile
    
    day_n_forward -= 1 ;; Remove 1 day to compensate that when the loop exits after finding a file, each counter is +1 from the actual day offset that produced the match.
    
    ;; Backward in time (for 2 days max - arbitrary choice)
    day_n_backward=1
    count_backward=0
    while (day_n_backward le day_limit) and (count_backward eq 0) do begin
      
      this_start_time_backward = anytim(anytim(start_time) - day_n_backward * day_in_s, /vms)
      this_end_time_backward   = anytim(anytim(end_time) - day_n_backward * day_in_s, /vms)

      date_path_backward = get_fid(this_start_time_backward,this_end_time_backward,/full,delim='/')

      path  = type_path + date_path_backward[0] +'/CAL'
      found_files_backward = sock_find(site,filter,path=path,count=count_backward)

      day_n_backward += 1

    endwhile
    
    day_n_backward -= 1 ;; Remove 1 day to compensate that when the loop exits after finding a file, each counter is +1 from the actual day offset that produced the match.
    
    if (count_forward eq 0) and (count_backward eq 0) then begin
      
      message, [" ", " ", "No STIX calibration file was found within 2 days before or after " +date_path[0]+". Please, download a calibration file manually.", " ", " "] 
      
    endif else begin
      
      diff_days = day_n_forward - day_n_backward 
      ;; if the difference is negative, it means that the closest file is the one found forward in time
      ;; if the difference is positive, it means that the closest file is the one found backward in time
      ;; if the difference is zero, we arbitrarily select the file found forward in time
      
      if diff_days le 0 then begin
        
        count = count_forward
        found_files = found_files_forward
        date_path = date_path_forward
        
      endif else begin
        
        count = count_backward
        found_files = found_files_backward
        date_path = date_path_backward
        
      endelse
      
    endelse
    
  endif

  message, [" ", " ", "Download STIX calibration file recorded at " +date_path[0], " ", " "], /continue

  len_path = STRLEN(site+path)

  if count eq 1 then begin

    selected_file = found_files[0]
    
    ;; Extract the calibration file start time from the filename. This is used below to verify that the ELUT onboard at the time of calibration matches the one used during the input time range.
    
    len_full_path = STRLEN(found_files[0])
    filename = STRMID(found_files[0], len_path+1, len_full_path)

    string_dates = STRMID(filename, STRLEN('solo_CAL_stix-cal-energy_'), STRLEN(filename)-STRLEN('solo_CAL_stix-cal-energy_')-9) ;; 9 is the number of characters in the string '_V02.fits'

    mid_part = STRPOS(string_dates, '-')

    this_start_time_file = STRMID(string_dates, 0, mid_part)
    year  = STRMID(this_start_time_file, 0, 4)
    month = STRMID(this_start_time_file, 4, 2)
    day   = STRMID(this_start_time_file, 6, 2)
    hour  = STRMID(this_start_time_file, 9, 2)
    min   = STRMID(this_start_time_file, 11, 2)
    sec   = STRMID(this_start_time_file, 13, 2)
    start_time_file = year + '-' + month + '-' + day + 'T' + hour + ':' + min + ':' + sec

  endif else begin

    message, [" ", " ", "More than 1 STIX calibration file found for this day ("+date_path[0]+"). Download the one with largest duration.", " ", " "], /continue
    
    start_time_file = []
    end_time_file = []

    ;; Extract file names
    for i = 0,n_elements(found_files)-1 do begin
      
      len_full_path = STRLEN(found_files[i])
      filename = STRMID(found_files[i], len_path+1, len_full_path)

      string_dates = STRMID(filename, STRLEN('solo_CAL_stix-cal-energy_'), STRLEN(filename)-STRLEN('solo_CAL_stix-cal-energy_')-9) ;; 9 is the number of characters in the string '_V02.fits'

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
    
    start_time_file = start_time_file[idx_file]

  endelse
  
  
  ;; Check that the ELUT onboard when the calibration file was recorded matches the ELUT onboard during the input time range; otherwise, raise an error.

  ;; Compare ELUT
  elut_file = stx_date2elut_file(start_time_file)
  elut_time_range = stx_date2elut_file(start_time)
  if elut_file ne elut_time_range then message, $
    [" ", " ", "The closest-in-time calibration file was acquired with an onboard ELUT that differs from the one used during the selected time range. Please, manually download a different calibration file.", " ", " "]

  sock_copy, selected_file, out_name, local_file=out_file, out_dir = out_dir, clobber=clobber
  
  filename = STRMID(out_file, STRLEN(out_dir), STRLEN(out_file)-STRLEN(out_dir))
  
  message, [" ", " ", "Downloaded file: " + filename, " ", " "], /continue
  
  return, out_file

  

end
