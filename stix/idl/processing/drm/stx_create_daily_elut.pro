;+
;
; NAME:
;
;   stx_create_daily_elut
;
; PURPOSE:
;
;   Returns the low and high energy bin edges from the ELUT closest in time to the flaring event.
;
; CALLING SEQUENCE:
;
;   daily_elut = stx_create_daily_elut(flare_time)
;
; INPUTS:
;
;   flare_time: string containing the start time of the flare. It must be compatible with the anytim function
;
; OUTPUTS:
;
;   Structure containing:
;   - ENERGY_BIN_LOW: float array of dimension 32 x 12 x 32 (number of energy bins x number of pixels x number of detectors) 
;                     containing the low energy edges of the daily ELUT.
;   - ENERGY_BIN_HIGH: float array of dimension 32 x 12 x 32 (number of energy bins x number of pixels x number of detectors) 
;                     containing the high energy edges of the daily ELUT.
;
; KEYWORDS:
;
;  gain_offset_version: string indicating which version of the daily gain/offset values is used. Options are:
;                   'original': values derived from O. Limousin's calibration runs
;                   'savgol': values smoothed by means of the Savitzky–Golay filter. The considered window size is 21 days.
;                   'median': values smoothed by means of a moving median filter. The considered window size is 21 days. This is the default value.
;
; HISTORY: 
;
; CONTACT:
;   paolo.massa@fhnw.ch
;-

function stx_create_daily_elut, flare_time, gain_offset_version=gain_offset_version

default, gain_offset_version, 'median'

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

idx_today = value_locate(calibration_date, anytim(flare_time))
today_gain = reform(daily_gain[*,*,idx_today])
today_offset = reform(daily_offset[*,*,idx_today])

;; Read ELUT channels
elut_filename = stx_date2elut_file(stx_time2any(flare_time))
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


return, {energy_bin_low: energy_bin_low, $
         energy_bin_high: energy_bin_high}

end