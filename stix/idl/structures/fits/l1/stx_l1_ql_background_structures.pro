;+
; :project:
;       STIX
;
; :name:
;       stx_l1_ql_background_structures
;
; :purpose:
;       Makes structures to store quick look l1 backgound TM and serialise to fits
;
; :categories:
;       structures, telemetry
;
; :params:
;    n_energies : in, type="int"
;             Number of energies channels
;
;    n_samples : in, type="int"
;             Number of variance measurements
;
; :returns:
;    A dictionary containing the energy and count structures
;
; :examples:
;    structures = stx_l1_ql_background_structures(5, 10)
;    control_struc = structures.energy
;    data_struc = structures.count
;
; :history:
;       02-May-2018 – SAM (TCD) init
;
;-
function stx_l1_ql_background_structures, n_energies, n_times

    energy = {CHANNEL: 0L, E_MIN: 0.0, E_MAX: 0.0}
    
    energies = replicate(energy, n_energies)

    count = {BACKGROUND: lonarr(n_energies), TRIGGERS: 0L, CHANNEL: lonarr(n_energies), TIME: 0.0d, $
        TIMEDEL: 0.0, LIVETIME: 1, ERROR: lonarr(n_energies)}
        
    counts = replicate(count, n_times)
    
    return, DICTIONARY('energy', energies, 'count', counts)
end