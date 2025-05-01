;+
;
; NAME:
;
;   stx_construct_pixel_data
;
; PURPOSE:
;
;   Read a STIX science L1 fits file (and optionally a STIX background L1 fits file) and fill in a 'pixel_data' structure
;
; CALLING SEQUENCE:
;
;   pixel_data = stx_construct_pixel_data(path_sci_file, time_range, energy_range)
;
; INPUTS:
;
;   path_sci_file: path of the STIX science L1 fits file
;
;   time_range: string array containing the start and the end of the time interval to consider
;
;   energy_range: two-element array containing the lower and upper edges of the energy interval to consider
;
; OUTPUTS:
;
;   'stx_pixel_data' structure containing:
;
;     - COUNTS: array 32x12 containing the number of counts recorded by the detector pixels
;               in the selected time and energy intervals
;     - COUNTS_ERROR: array 32x12 containing the errors (statistics + compression) associated with the number of counts
;                     recorded by the detector pixels
;     - LIVE_TIME: 32-element array containing the live time of each detector in the considered time interval
;     - TIME_RANGE: two-element 'stx_time' array containing the lower and upper edge of the selected time interval
;                   (the time bins containing the start and the end time provided as input are included in the
;                   selected interval)
;     - ENERGY_RANGE: array containing the lower and the upper edge of the selected energy interval
;                     (the energy bins containing the lower and the upper edges provided as input are included in the interval)
;     - TOT_COUNTS_BKG: double, estimate of the total number of background counts recorded during the flaring event
;     - XY_FLARE: two-element array containing the X and Y coordinates of the estimated flare location (STIX coordinate frame, arcsec).
;                   If 'xy_flare' is not passed, it is filled with NaN values
;     - RCR: rate control regime status in the selcted time interval. If the RCR changes in that interval, an error is thrown
;     - PIXEL_MASKS: matrix containing info on the pixels used in the selected time interval
;     - DETECTOR_MASKS: matrix containing information on the detectors used in the selected time interval
;
; KEYWORDS:
;
;   path_bkg_file: path of a background L1 fits file. If provided, the fields 'COUNTS_BKG', 'COUNTS_ERROR_BKG' and 'LIVE_TIME_BKG'
;                  of the pixel_data structure are filled with the values read from the background measurement file
;
;   elut_corr: if set, a correction based on a ELUT table is applied to the measured counts
;
;   xy_flare: two-element array containing the X and Y coordinates of the estimated flare location
;             (STIX coordinate frame, arcsec). If passed, the grid transmission correction is computed
;
;   subc_index: array containing the indices of the selected imaging detectors. Used only for plotting the lightcurve by means of
;             'stx_plot_selected_time_range'. Default, indices of
;             the detectors from 10 to 3
;
;   sumcase: string indicating which pixels are considered. Used only for plotting the lightcurve by means of
;            'stx_plot_selected_time_range'. Refer to the header of 'stx_sum_pixel_data' for more details.
;             Default, 'TOP+BOT'
;
;   silent: if set, plots are not displayed
;
;   no_small: if set, Moire patterns measured by small pixels are not plotted with 'stx_plot_moire_pattern'
;
;   no_rcr_check: if set, control on RCR change during the selected time interval is not performed. Default, 0
;
;   shift_duration: if set, shift all time bins by 1 to account for FSW time input discrepancy prior to 09-Dec-2021.
;                   N.B. WILL ONLY WORK WITH FULL TIME RESOLUTION DATA WHICH IS OFTEN NOT THE CASE FOR PIXEL DATA.
;
;   gain_offset_version: string indicating which version of the daily gain/offset values is used. Options are:
;                   'original': values derived from O. Limousin's calibration runs
;                   'savgol': values smoothed by means of the Savitzky–Golay filter. The considered window size is 21 days.
;                   'median': values smoothed by means of a moving median filter. The considered window size is 21 days. This is the default value.
;
; HISTORY: July 2022, Massa P., created
;          September 2022, Massa P., added 'shift_duration' keyword
;          May 2023, Massa P., do not call 'stx_plot_selected_time_range' if the science fits file contains a single time bin
;          October 2023, Massa P., fixed bug in the selection of the energy bin indices
;          November 2023, Massa P., use simplified version of the subcollimator transmission (temporary solution)
;          April, 2025, Massa P., new ELUT correction based on daily gain and offset values. 
;              Bkg subtraction is applied before ELUT correction is computed.
;
; CONTACT:
;   paolo.massa@fhnw.ch
;-

function stx_construct_pixel_data, path_sci_file, time_range, energy_range, elut_corr=elut_corr, $
  path_bkg_file=path_bkg_file, xy_flare=xy_flare, subc_index=subc_index, $
  sumcase=sumcase, silent=silent, no_small=no_small, no_rcr_check=no_rcr_check, $
  shift_duration=shift_duration, gain_offset_version=gain_offset_version,  _extra=extra

  default, elut_corr, 1
  default, silent, 0
  default, subc_index, stx_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c',$
    '6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c'])

  default, sumcase, "TOP+BOT"
  default, no_rcr_check, 0
  default, gain_offset_version, 'median'

  if anytim(time_range[0]) gt anytim(time_range[1]) then message, "Start time is greater than end time"
  if energy_range[0] gt energy_range[1] then message, "Energy range lower edge is greater than the higher edge"


  stx_read_pixel_data_fits_file, path_sci_file, data_str = data, t_axis = t_axis, e_axis = e_axis, alpha = alpha, shift_duration=shift_duration, _extra=extra

  if keyword_set(path_bkg_file) then begin

    stx_read_pixel_data_fits_file, path_bkg_file, data_str = data_bkg, t_axis = t_axis_bkg, e_axis = e_axis_bkg, _extra=extra

    if n_elements(t_axis_bkg.DURATION) gt 1 then message, 'The chosen file does not contain a background measurement'

  endif


  ;;************** Select time indices

  time_start = stx_time2any(t_axis.TIME_START)
  time_end   = stx_time2any(t_axis.TIME_END)
  if (anytim(time_range[0]) lt min(time_start)) then $
    message, 'The selected start time is outside the time interval of the data file (' + $
    anytim(min(time_start), /vms) + ' - ' + anytim(max(time_end), /vms) + ')'

  if (anytim(time_range[1]) gt max(time_end)) then $
    message, 'The selected end time is outside the time interval of the data file (' + $
    anytim(min(time_start), /vms) + ' - ' + anytim(max(time_end), /vms) + ')'

  time_ind_min = where(time_start le anytim(time_range[0]))
  time_ind_min = time_ind_min[-1]
  time_ind_max = where(time_end ge anytim(time_range[1]))
  time_ind_max = time_ind_max[0]

  time_ind        = [time_ind_min:time_ind_max]
  this_time_range = [t_axis.TIME_START[time_ind_min], t_axis.TIME_END[time_ind_max]]

  ;;************** Select energy indices

  ;; Select indices of the energy bins (among the 32) that are actually present in the pixel data science file
  energy_bin_mask = data.energy_bin_mask
  energy_bin_idx  = where(energy_bin_mask eq 1)
  
  energy_low  = e_axis.LOW
  energy_high = e_axis.HIGH
  
  if (energy_range[0] lt min(energy_low)) then $
    message, 'The lower edge of the selected energy interval is outside the science energy interval (' + $
    num2str(fix(min(energy_low))) + ' - ' + num2str(fix(max(energy_high))) + ' keV)'

  if (energy_range[1] gt max(energy_high)) then $
    message, 'The upper edge of the selected energy interval is outside the science energy interval (' + $
    num2str(fix(min(energy_low))) + ' - ' + num2str(fix(max(energy_high))) + ' keV)'

  energy_min = min(energy_low)
  energy_max = max(energy_high)

  if keyword_set(path_bkg_file) then begin

    ;; Select indices of the energy bins (among the 32) that are actually present in the pixel data bkg file
    energy_bin_mask_bkg = data_bkg.energy_bin_mask
    energy_bin_idx_bkg = where(energy_bin_mask_bkg eq 1)

    energy_low_bkg  = e_axis_bkg.LOW
    energy_high_bkg = e_axis_bkg.HIGH

    if (energy_range[0] lt min(energy_low_bkg)) then $
      message, 'The lower edge of the selected energy interval is outside the background energy interval (' + $
      num2str(fix(min(energy_low_bkg))) + ' - ' + num2str(fix(max(energy_high_bkg))) + ' keV)'

    if (energy_range[1] gt max(energy_high_bkg)) then $
      message, 'The upper edge of the selected energy interval is outside the background energy interval (' + $
      num2str(fix(min(energy_low_bkg))) + ' - ' + num2str(fix(max(energy_high_bkg))) + ' keV)'

    ;; Extract energy range in common between science and background file
    energy_min = max([energy_min,min(energy_low_bkg)])
    energy_max = min([energy_max,max(energy_high_bkg)])
    idx_energy_bkg = where((energy_low_bkg ge energy_min) and (energy_high_bkg le energy_max))
    
    energy_bin_idx_bkg = energy_bin_idx_bkg[idx_energy_bkg]
    energy_low_bkg = energy_low_bkg[idx_energy_bkg]
    energy_high_bkg = energy_high_bkg[idx_energy_bkg]
    
    energy_ind_min_bkg = where(energy_low_bkg le energy_range[0])
    energy_ind_min_bkg = energy_ind_min_bkg[-1]
    energy_ind_max_bkg = where(energy_high_bkg ge energy_range[1])
    energy_ind_max_bkg = energy_ind_max_bkg[0]

    energy_ind_bkg        = [energy_ind_min_bkg:energy_ind_max_bkg]
    if n_elements(energy_ind_bkg) eq 1 then energy_ind_bkg = energy_ind_bkg[0]

    this_energy_range_bkg = [energy_low_bkg[energy_ind_min_bkg], energy_high_bkg[energy_ind_max_bkg]]

  endif
  
  ;; Extract energy range in common between science and background file
  idx_energy = where((energy_low ge energy_min) and (energy_high le energy_max))
  
  energy_bin_idx = energy_bin_idx[idx_energy]
  energy_low = energy_low[idx_energy]
  energy_high = energy_high[idx_energy]
  
  energy_ind_min = where(energy_low le energy_range[0])
  energy_ind_min = energy_ind_min[-1]
  energy_ind_max = where(energy_high ge energy_range[1])
  energy_ind_max = energy_ind_max[0]

  energy_ind        = [energy_ind_min:energy_ind_max]
  if n_elements(energy_ind) eq 1 then energy_ind = energy_ind[0]

  this_energy_range = [energy_low[energy_ind_min], energy_high[energy_ind_max]]

  ;;************** Compute livetime

  triggergram        = stx_triggergram(data.TRIGGERS, data.TRIGGERS_ERR, t_axis)
  livetime_fraction_data = stx_livetime_fraction(triggergram)
  livetime_fraction = livetime_fraction_data.livetime_fraction
  livetime_fraction_err = livetime_fraction_data.livetime_fraction_err
  duration_time_bins = t_axis.DURATION
  duration_time_bins = transpose(cmreplicate(duration_time_bins, 32))
  
  live_time_bins     = livetime_fraction * duration_time_bins
  live_time_bins_err = livetime_fraction_err * duration_time_bins
  
  live_time = n_elements(time_ind) eq 1? reform(live_time_bins[*,time_ind]) : total(live_time_bins[*,time_ind],2)
  live_time_error = n_elements(time_ind) eq 1? reform(live_time_bins_err[*,time_ind]) : $
    sqrt(total(live_time_bins_err[*,time_ind]^2.,2))
  
  if keyword_set(path_bkg_file) then begin
  
    triggergram_bkg        = stx_triggergram(data_bkg.TRIGGERS, data_bkg.TRIGGERS_ERR, t_axis_bkg)
    livetime_fraction_bkg_data = stx_livetime_fraction(triggergram_bkg)
    livetime_fraction_bkg = livetime_fraction_bkg_data.livetime_fraction
    livetime_fraction_bkg_err = livetime_fraction_bkg_data.livetime_fraction_err
    duration_time_bins_bkg = t_axis_bkg.DURATION
    live_time_bkg          = duration_time_bins_bkg[0]*livetime_fraction_bkg
    live_time_error_bkg    = duration_time_bins_bkg[0]*livetime_fraction_bkg_err
  
  endif else begin
  
    live_time_bkg = dblarr(32) + 1.
    live_time_error_bkg = dblarr(32)
    
  endelse

  ;;************** Define count matrix and bkg count matrix

  ;; WE ASSUME THAT THE SCIENCE DATA FILE AND THE BKG DATA FILE CONTAIN MORE THAN 1 ENERGY BIN
  ;; (I.E., THE NUMBER OF ELEMENTS IN energy_bin_idx AND IN energy_bin_idx_bkg IS ASSUMED TO BE
  ;; LARGER THAN 1)

  if (n_elements(energy_bin_idx) eq 1) then $
    message, 'It is not possible to analyze this science file as it contains a single energy bin. Please, contact the STIX team for further details.'
  if keyword_set(path_bkg_file) and (n_elements(energy_bin_idx_bkg) eq 1) then $
    message, 'It is not possible to utilize this bkg file since it contains a single energy bin. Please, either utilize a different bkg file or contact the STIX team.'

  ;; Dimensions: [energy,pixel,detector,time]
  counts       = data.COUNTS
  counts_error = data.COUNTS_ERR
  
  ;; Consider only selected energy bins
  counts       = counts[energy_bin_idx,*,*,*]
  counts_error = counts_error[energy_bin_idx,*,*,*]

  if keyword_set(path_bkg_file) then begin
    
    ;; Check if science and background files are reconrded with the same ELUT
    
    elut_filename = stx_date2elut_file(stx_time2any(this_time_range[0]))
    stx_read_elut, elut_gain, elut_offset, adc4096_str, elut_filename = elut_filename, ekev_actual = ekev_actual_elut
    
    elut_filename_bkg = stx_date2elut_file(stx_time2any(t_axis_bkg.TIME_START))
    stx_read_elut, ekev_actual = ekev_actual_bkg, elut_filename = elut_filename_bkg
    
    ;; Compare ELUT tables
    elut_comp = STRCMP(elut_filename, elut_filename_bkg)
    
    if not elut_comp then $
      message, 'The background file must be recorded when the same ELUT as the science file was uploaded. Please choose a different background file that is closer in time to the science file.'

    counts_bkg       = data_bkg.COUNTS
    counts_error_bkg = data_bkg.COUNTS_ERR
    counts_bkg       = counts_bkg[energy_bin_idx_bkg,*,*]
    counts_error_bkg = counts_error_bkg[energy_bin_idx_bkg,*,*]

  endif else begin
    
    counts_bkg       = dblarr(size(counts, /dim))
    counts_error_bkg = dblarr(size(counts, /dim))
    
  endelse
  
  ;;************** Plot lightcurve (if ~silent)

  if ~silent and (n_elements(size(live_time_bins, /dim)) gt 1) then begin

    if keyword_set(path_bkg_file) then begin
      stx_plot_selected_time_range, stx_time2any(t_axis.MEAN), energy_ind, time_ind, counts, live_time_bins, subc_index, $
        sumcase, this_energy_range, this_time_range, counts_bkg=counts_bkg, $
        live_time_bkg=live_time_bkg, energy_ind_bkg=energy_ind_bkg
    endif else begin
      stx_plot_selected_time_range, stx_time2any(t_axis.MEAN), energy_ind, time_ind, counts, live_time_bins, subc_index, $
        sumcase, this_energy_range, this_time_range
    endelse

  endif

  ;;************** Sum counts (and related errors) in time

  if n_elements(time_ind) eq 1 then begin

    counts       = reform(counts[*,*,*,time_ind])
    counts_error = reform(counts_error[*,*,*,time_ind])

  endif else begin

    counts       = total(counts[*,*,*,time_ind],4)
    counts_error = sqrt(total(counts_error[*,*,*,time_ind]^2.,4))

  endelse
  
  ;;************** Background subtraction
  
  live_time_rep = transpose(cmreplicate(live_time, [n_elements(energy_bin_idx),12]), [1,2,0])
  live_time_error_rep = transpose(cmreplicate(live_time_error, [n_elements(energy_bin_idx),12]), [1,2,0])
  live_time_bkg_rep = keyword_set(path_bkg_file)? transpose(cmreplicate(live_time_bkg, [n_elements(energy_bin_idx_bkg),12]), [1,2,0]) : fltarr(size(counts_bkg, /dim)) + 1.
  live_time_error_bkg_rep = keyword_set(path_bkg_file)? transpose(cmreplicate(live_time_error_bkg, [n_elements(energy_bin_idx_bkg),12]), [1,2,0]) : fltarr(size(counts_bkg, /dim))
  
  counts_bkg_estimate = f_div( live_time_rep * counts_bkg, live_time_bkg_rep )
  error_numerator = abs(live_time_rep * counts_bkg) * sqrt( f_div(live_time_error_rep,live_time_rep)^2. + $
    f_div(counts_error_bkg,counts_bkg)^2.)
  counts_bkg_estimate_error = abs(counts_bkg_estimate) * sqrt( f_div(error_numerator,live_time_rep * counts_bkg)^2. + $
    f_div(live_time_error_bkg_rep,live_time_bkg_rep)^2.)
  
  if elut_corr and ~silent then begin

    ;; To be used for plot of the bkg subtracted spectrum
    spectrum_with_bkg = total(total(counts[*,0:7,subc_index], 3), 2) / (energy_high - energy_low)
    spectrum_bkg = keyword_set(path_bkg_file)? total(total(counts_bkg_estimate[*,0:7,subc_index], 3), 2) / (energy_high_bkg - energy_low_bkg) : total(total(counts_bkg[*,0:7,subc_index], 3), 2)
    
  endif
  
  ;; Determine indices of the selected pixels. It will be used for estimating total number of BKG counts and for estimating the amount of ELUT correction
  case sumcase of

    'TOP':     begin
      pixel_ind = [0]
    end

    'BOT':     begin
      pixel_ind = [1]
    end

    'TOP+BOT': begin
      pixel_ind = [0,1]
    end

    'ALL': begin
      pixel_ind = [0,1,2]
    end

    'SMALL': begin
      pixel_ind = [2]
    end
  end
  
  ;; Compute total number of counts (to be used for comparison with bkg counts)
  counts_reshaped = reform(counts,n_elements(energy_bin_idx), 4, 3, 32)
  tot_counts = total(counts_reshaped[energy_ind,*,pixel_ind,subc_index])
  
  ;; Apply BKG subtraction
  counts = counts - counts_bkg_estimate
  counts_error = sqrt(counts_error^2. +  counts_bkg_estimate_error^2.)  
  
  if keyword_set(path_bkg_file) then begin
    
    ;; Compute total number of background counts. Select only imaging detectors
    counts_bkg = reform(counts_bkg_estimate, n_elements(energy_bin_idx_bkg), 4, 3, 32)
    tot_counts_bkg = total(counts_bkg[energy_ind_bkg,*,pixel_ind,subc_index])
  
  endif
  
  ;;************** Print total number of counts

  if ~silent then begin

    print
    print
    print,'***********************************************************************'
    print,'Total number of counts in image:  '+strtrim(tot_counts)
    print,'Background counts:                '+strtrim(tot_counts_bkg)
    print,'Counts above background:          '+strtrim(tot_counts-tot_counts_bkg)
    print,'Total to background:              '+strtrim(tot_counts/tot_counts_bkg)
    print,'***********************************************************************'
    print
    print

  endif
  
  ;;************** Sum counts (and related errors) in energy

  ;; Compute elut correction (if 'elut_corr' is set) - Correct just the first and the last energy bins (flat spectrum is assumed)
  if elut_corr then begin
    
    ;;**************** Compute spectral index to be used for ELUT correction
    
    spectrum = total(total(counts[*,0:7,subc_index], 3), 2) / (energy_high - energy_low)
    index_data = stx_estimate_spectral_index(energy_low, energy_high, spectrum)
    
    sp_index = index_data.index_final
    idx_peak = index_data.idx_peak
    
    if ~silent then begin
    
      charsize = 1.8
      
      loadct,5
      device, Window_State=win_state
      if not win_state[10] then window,10,xsize=1000,ysize=500
      wset,10
      cleanplot
      
      !p.multi = [0,2,1]
      
      energy_axis = (energy_low+energy_high)/2.
      
      plot, energy_axis, spectrum_with_bkg, psym=10, /xst, /yst, /ylog, /xlog, yrange=[max([1.,min(spectrum)]),max(spectrum)*10.],charsize=charsize, $
        title='STIX spectrum',xtitle='Energy [keV]', ytitle = 'STIX spectrum [counts s!U-1!n keV!U-1!n]'
      oplot, [this_energy_range[0],this_energy_range[0]], [max([1.,min(spectrum)]),max(spectrum)*10.], linestyle=1
      oplot, [this_energy_range[1],this_energy_range[1]], [max([1.,min(spectrum)]),max(spectrum)*10.], linestyle=1
      oplot, energy_axis, spectrum_bkg, psym=10, linestyle=2
      oplot, energy_axis, spectrum, psym=10, color=122
      leg_text = ['Observed spectrum', 'Background', 'BKG-subtracted']
      leg_color = [255,255,122]
      leg_style = [0, 2, 0]
      ssw_legend, leg_text, color=leg_color, linest=leg_style, box=0, charsize=1.5, thick=1.5, /right
      
      
      plot,energy_axis, sp_index, psym=10, /xst, /yst, /xlog, charsize=charsize, $
        title='Estimate of the spectral index', xtitle='Energy [keV]', ytitle = 'Spectral index'
      oplot, [this_energy_range[0],this_energy_range[0]], [-20,20], linestyle=1
      oplot, [this_energy_range[1],this_energy_range[1]], [-20,20], linestyle=1
      oplot, [4,150], [0,0], linestyle=2
      
    endif
    
    ;; Read daily gain and offset
    calibration_data = mrdfits(concat_dir( concat_dir('SSW_STIX','dbase'),'detector') + get_delim() + 'daily_gain_offset.fits',1)
    calibration_date = anytim(calibration_data.DATE)
    
    case gain_offset_version of
    
      'original': begin
        
                  daily_gain = calibration_data.GAIN
                  daily_offset = calibration_data.OFFSET
                  
                  end

      'savgol': begin
                
                daily_gain = calibration_data.GAIN_SAVGOL
                daily_offset = calibration_data.OFFSET_SAVGOL
                
                end

      'median': begin
        
                daily_gain = calibration_data.GAIN_MEDIAN
                daily_offset = calibration_data.OFFSET_MEDIAN
                
                end
                
        else: begin
              message, 'Keyword gain_offset_version must be set equal to original, savgol, or median.'
              end
    
    end
    
    
    if anytim(time_range[0]) le calibration_date[0] then begin
      
      print
      print
      print, 'WARNING: daily calibration data are available only from ' + anytim(calibration_date[0], /vms) + '. Results can be inaccurate!'
      print
      print
      
    endif
    
    idx_today = value_locate(calibration_date, anytim(time_range[0]))
    today_gain = reform(daily_gain[*,*,idx_today])
    today_offset = reform(daily_offset[*,*,idx_today])

    ;; Read ELUT channels
    elut_filename = stx_date2elut_file(stx_time2any(this_time_range[0]))
    stx_read_elut, elut_gain, elut_offset, adc4096_str, elut_filename = elut_filename, ekev_actual = ekev_actual_elut
    
    
    elut_channels = adc4096_str.ADC4096
    ekev_actual = fltarr(31,12,32)
    
    for i=0,30 do begin
      ;; Energy edges in keV which correspond to the onboard binning edges in native channel unit (~0.1 keV)
      ;; ekev_actual is obtained by applying the daily gain and offset derived from calibration fitting to the onboard ELUT channel numbers
      ekev_actual[i,*,*] =  (reform(elut_channels[i,*,*]) - today_offset) / today_gain 
      
    endfor
    
    ;; We assume that the first energy bin starts from 0. The right end of the last energy bin is set to NaN
    ;; We assume that the energy interval is always partitioned into 32 energy bins
    energy_bin_low           = fltarr(32,12,32)
    energy_bin_low[1:31,*,*] = ekev_actual

    energy_bin_high           = fltarr(32,12,32)
    energy_bin_high[0:30,*,*] = ekev_actual
    energy_bin_high[31,*,*]   = !VALUES.F_NaN

    energy_bin_low  = energy_bin_low[energy_bin_idx,*,*]
    energy_bin_high = energy_bin_high[energy_bin_idx,*,*]
    

    ;; We assume the spectrum has a powerlaw distribution E^-sp_index at any energy bin.
    ;
    ;  The number of counts in the energy range [a,b] is 
    ;  
    ;  C = \int_a^b E^-sp_index dE = (b^(-sp_index+1) - a^(-sp_index+1)) / (-sp_index+1)
    ;
    ;  We distinguish 3 cases:
    ;
    ; 1. The considered energy range consists of a single energy bin which is not at the peak 
    ;    (i.e., different from 6-7 keV; example: 9-10 keV). The correction factor is
    ;    
    ;    corr = (10^(-sp_index+1) - 9^(-sp_index+1)) / (10.03^(-sp_index+1) - 8.97^(-sp_index+1)), 
    ;    
    ;    where 10.03 and 8.97 are derived from daily ELUT
    ;
    ; 2. The considered energy range consists of a single energy bin which is at the peak of the spectrum 
    ;    (i.e., 6-7 keV when attenuator is not inserted). We consider the middle point of the energy range (e.g., 6.5 keV) 
    ;    and we assume that the spectrum follows a powerlaw distribution below and above the middle point, viz.
    ;    
    ;             A E^-alpha if E < 6.5 keV
    ;     F(E) =  
    ;             B E^-beta if E >= 6.5 keV
    ;    
    ;    where alpha < 0 and beta > 0. Assuming that the spectrum is continuous at 6.5 keV, 
    ;    we obtain that A and B are related to the following equation:
    ;    
    ;    A = B 6.5^(alpha - beta)
    ;    
    ;    Assuming that the ELUT edges of the considered energy bin are 5.98 and 7.05 keV, we obtain
    ;    
    ;    corr = (\int_6^7 F(E) dE) / (\int_5.98^7.05 F(E) dE) ,
    ;    
    ;    which leads to
    ;    
    ;    corr = (6.5^(alpha - beta) * (6.5^(-alpha+1) - 6^(-alpha+1)) / (-alpha + 1) ) + (7^(-beta+1) - 6.5^(-beta+1)) / (-beta+1) ) / 
    ;           (6.5^(alpha - beta) * (6.5^(-alpha+1) - 5.98^(-alpha+1)) / (-alpha + 1) ) + (7.05^(-beta+1) - 6.5^(-beta+1)) / (-beta+1) )
    ;    
    ; 3. The considered energy range consists of multiple energy bins (e.g, 5-10 keV). Therefore, we apply a correction factor 
    ;    only to the lowest and to the highest bin (corr_low and corr_high). Assuming that the energy edges of these bins contained 
    ;    in the ELUT table are 5.01, 5.98 keV and 8.97, 10.02 keV, we have
    ;    
    ;    corr_low = (5.98^(-alpha+1) - 5^(-alpha+1)) / (5.98^(-alpha+1) - 5.01^(-alpha+1))
    ;    
    ;    and
    ;    
    ;    corr_high = (10^(-beta+1) - 8.97^(-beta+1)) / (10.02^(-beta+1) - 8.97^(-beta+1))
    ;    
    ;    where alpha and beta are the powerlaw indices in the first and the last energy bins, respectively.
  
    if n_elements(energy_ind) eq 1 then begin
      
      if energy_ind eq idx_peak then begin
        
        ;; energy_bin_idx contains at least two indices (see control above)
        case energy_ind of
          
          0: begin
            
            sp_index_low = sp_index[energy_ind+1]
            sp_index_high = sp_index[energy_ind+1]
            
            end
           
          n_elements(energy_bin_idx)-1: begin
            
            sp_index_low = sp_index[energy_ind-1]
            sp_index_high = sp_index[energy_ind-1]
            
            end
            
          else: begin
            
            sp_index_low = sp_index[energy_ind-1]
            sp_index_high = sp_index[energy_ind+1]
            
            end
          
        endcase
        
        energy_bin_mid = (energy_high[energy_ind] + energy_low[energy_ind]) / 2.
        
        corr_factor_low_ELUT = (energy_bin_mid^(-sp_index_low+1.) - reform(energy_bin_low[energy_ind,*,*])^(-sp_index_low+1.)) / $
                               (-sp_index_low+1.)
        corr_factor_high_ELUT = (reform(energy_bin_high[energy_ind,*,*])^(-sp_index_high+1.) - energy_bin_mid^(-sp_index_high+1.)) / $
                               (-sp_index_high+1.)
        
        corr_factor_low_SCI = (energy_bin_mid^(-sp_index_low+1.) - energy_low[energy_ind]^(-sp_index_low+1.)) / (-sp_index_low+1.)
        corr_factor_high_SCI = (energy_high[energy_ind]^(-sp_index_high+1.) - energy_bin_mid^(-sp_index_high+1.)) / (-sp_index_high+1.)
  
        norm_factor = energy_bin_mid^(-sp_index_high + sp_index_low)
        energy_corr_factor = elut_corr ? (norm_factor * corr_factor_low_SCI + corr_factor_high_SCI) / $
                                         (norm_factor * corr_factor_low_ELUT + corr_factor_high_ELUT) : dblarr(12,32)+1
        
      endif else begin
  
        energy_corr_factor = elut_corr ? $
          (energy_high[energy_ind]^(-sp_index[energy_ind]+1.) - energy_low[energy_ind]^(-sp_index[energy_ind]+1.)) / $
          (reform(energy_bin_high[energy_ind,*,*])^(-sp_index[energy_ind]+1.) - reform(energy_bin_low[energy_ind,*,*])^(-sp_index[energy_ind]+1.)) : $
          dblarr(12,32)+1
  
      endelse
      
      ;; Compute total number of counts BEFORE ELUT correction is applied. It is used for estimating the amount of the ELUT correction
      counts_reshaped = reform(counts, n_elements(energy_bin_idx), 4, 3, 32)
      counts_no_elut = reform(total(counts_reshaped[energy_ind,*,pixel_ind,subc_index], 1))
      
      counts     = reform(counts[energy_ind,*,*]) * energy_corr_factor
      counts_error = reform(counts_error[energy_ind,*,*]) * energy_corr_factor
      
      ;; Compute total number of counts AFTER ELUT correction is applied. It is used for estimating the amount of the ELUT correction
      counts_reshaped = reform(counts, 4, 3, 32)
      counts_elut = reform(counts_reshaped[*,pixel_ind,subc_index])
  
    endif else begin
    
      energy_corr_factor_low = elut_corr ? $
        (reform(energy_bin_high[energy_ind[0],*,*])^(-sp_index[energy_ind[0]]+1.) - energy_low[energy_ind[0]]^(-sp_index[energy_ind[0]]+1.)) / $
        (reform(energy_bin_high[energy_ind[0],*,*])^(-sp_index[energy_ind[0]]+1.) - reform(energy_bin_low[energy_ind[0],*,*])^(-sp_index[energy_ind[0]]+1.)) : $
        dblarr(12,32)+1
  
      energy_corr_factor_high = elut_corr ? $
        (energy_high[energy_ind[-1]]^(-sp_index[energy_ind[-1]]+1.) - reform(energy_bin_low[energy_ind[-1],*,*])^(-sp_index[energy_ind[-1]]+1.)) / $
        (reform(energy_bin_high[energy_ind[-1],*,*])^(-sp_index[energy_ind[-1]]+1.) - reform(energy_bin_low[energy_ind[-1],*,*])^(-sp_index[energy_ind[-1]]+1.)) : $
        dblarr(12,32)+1
  
      ;; Compute total number of counts BEFORE ELUT correction is applied. It is used for estimating the amount of the ELUT correction
      counts_reshaped = reform(counts, n_elements(energy_bin_idx), 4, 3, 32)
      counts_no_elut = reform(total(counts_reshaped[energy_ind,*,pixel_ind,subc_index], 1))
      
      counts[energy_ind[0],*,*]        *= energy_corr_factor_low
      counts[energy_ind[-1],*,*]       *= energy_corr_factor_high
      counts_error[energy_ind[0],*,*]  *= energy_corr_factor_low
      counts_error[energy_ind[-1],*,*] *= energy_corr_factor_high
  
      counts       = total(counts[energy_ind,*,*],1)
      counts_error = sqrt(total(counts_error[energy_ind,*,*]^2.,1))
      
      ;; Compute total number of counts AFTER ELUT correction is applied. It is used for estimating the amount of the ELUT correction
      counts_reshaped = reform(counts, 4, 3, 32)
      counts_elut = reform(counts_reshaped[*,pixel_ind,subc_index])
  
    endelse
  
    if ~silent then begin

       ;; Print minimun and maximum ELUT correction in the different pixels
      elut_corr_perc = f_div(abs(counts_no_elut -  counts_elut),counts_no_elut) * 100 ;; in percentage
      counts_error_reshaped = reform(counts_error, 4, 3, 32)
      counts_error_elut=reform(counts_error_reshaped[*,pixel_ind,subc_index])
      rel_error=abs(counts_error_elut/counts_elut)*100.
      print
      print
      print,'***********************************************************************'
      print,'Min/Max ELUT correction in the different pixels: '+num2str(min(elut_corr_perc), format='(f7.2)')+'% - '+num2str(max(elut_corr_perc), format='(f7.2)')+'% '
      print,'Averaged ELUT correction in a pixel:             '+num2str(average(elut_corr_perc), format='(f7.2)')+'%'
      print,'standard deviation of ELUT correction:           '+num2str(stdev(elut_corr_perc), format='(f7.2)')+'%'
      print,'Averaged error in counts in a pixel:             '+num2str(average(rel_error), format='(f7.2)')+'%'
      print,'***********************************************************************'
      print
      print
      
    endif
  endif else begin
    
    counts       = total(counts[energy_ind,*,*],1)
    counts_error = sqrt(total(counts_error[energy_ind,*,*]^2.,1))
    
  endelse
  ;;************** Correction for grid internal shadowing

  if keyword_set(xy_flare) then begin

    ;; Use simplified version of the grid transmission (temporary solution)
    subc_transmission     = stx_subc_transmission(xy_flare, /simple_transm)
    subc_transmission_bkg = stx_subc_transmission([0.,0.], /simple_transm)
    for i=0,31 do begin

      counts[*,i]       = counts[*,i]/subc_transmission[i]*0.25
      counts_error[*,i] = counts_error[*,i]/subc_transmission[i]*0.25

    endfor

  endif

  ;;**************  RCR

  rcr = data.rcr[time_ind]
  find_changes, rcr, index, state, count=count
  if (count gt 1) and (~no_rcr_check) then message, "RCR status changed in the selected time interval"

  ;;************** Pixel masks and detector masks

  if alpha then begin

    pixel_masks    = data.PIXEL_MASKS[*,*,time_ind]

  endif else begin

    pixel_masks    = data.PIXEL_MASKS[*,time_ind]

  endelse
  detector_masks = data.DETECTOR_MASKS[*,time_ind]

  diff_pixel_masks    = fltarr(n_elements(time_ind))
  diff_detector_masks = fltarr(n_elements(time_ind))

  for i=0,n_elements(time_ind)-2 do begin

    if alpha then begin

      diff_pixel_masks[i] = total(pixel_masks[*,*,i]-pixel_masks[*,*,i+1])

    endif else begin

      diff_pixel_masks[i] = total(pixel_masks[*,i]-pixel_masks[*,i+1])

    endelse

  endfor

  if total(diff_pixel_masks) gt 0. then message, "Pixel masks changed in the selected time interval"
  if total(diff_detector_masks) gt 0. then message, "Detector masks changed in the selected time interval"

  ;;************** Fill in pixel data structure

  pixel_data = stx_pixel_data()

  pixel_data.TIME_RANGE       = this_time_range
  pixel_data.ENERGY_RANGE     = this_energy_range
  pixel_data.LIVE_TIME        = live_time
  pixel_data.LIVE_TIME_ERROR  = live_time_error
  pixel_data.COUNTS           = transpose(counts)
  pixel_data.COUNTS_ERROR     = transpose(counts_error)
  pixel_data.TOT_COUNTS       = tot_counts
  if keyword_set(path_bkg_file) then pixel_data.TOT_COUNTS_BKG   = tot_counts_bkg

  if ~keyword_set(xy_flare) then begin
    pixel_data.XY_FLARE = [!VALUES.F_NaN,!VALUES.F_NaN]
  endif else begin
    pixel_data.XY_FLARE = xy_flare
  endelse

  pixel_data.RCR            = rcr[0]

  if alpha then begin
    pixels_used = where( total((data.pixel_masks)[*,*,0],1) eq 1, npix)
    pixel_data.PIXEL_MASKS[pixels_used] = 1b
  endif else begin
    pixel_data.PIXEL_MASKS    = pixel_masks[*,0]
  endelse

  pixel_data.DETECTOR_MASKS = detector_masks[*,0]

  if ~silent then stx_plot_moire_pattern, pixel_data, no_small=no_small

  return, pixel_data

end

