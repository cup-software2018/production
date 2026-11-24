# The integration windows of a run, from the DLY in its run config.  Source it,
# then call
#
#   fadc_windows <run config>    exports FADC_BASE_HI, FADC_INT_LO, FADC_INT_HI
#   iadc_windows <run config>    exports IADC_BASE_HI, IADC_INT_LO, IADC_INT_HI
#
# which prod_fadc.yaml and prod_iadc.yaml read as ${...}.  Each sets WINDOWS_NOTE
# to a one-line summary, or prints why and returns 1.  The run config is
# <CONF>/<run:06d>.yml.  Used by prod_wcmd_subrun.sh, prod_psmd_subrun.sh and
# autoprod, so the windows are decided here only.

# FADC, the same on every board.  At 2 ns per sample the light starts at
# E = 232 + (DLY - 1100) / 2 (main and door, runs 263, 301, 308, 342, to the
# sample); the baseline is [0, E - 32) and the integration [E, E + 184).
# E + 184 stays in the 1008-sample record up to DLY 2284.  The first DLY of the
# FADCT block.
fadc_windows() {
  local cfg=$1 dly edge
  [ -f "$cfg" ] || { echo "no run config $cfg for the FADC DLY"; return 1; }
  dly=$(awk '/^[^ -]/ { in_b = ($1 == "FADCT:"); next } in_b' "$cfg" |
        grep -oE '\bDLY: *[0-9]+' | grep -oE '[0-9]+' | head -1)
  [ -n "$dly" ] || { echo "no FADCT DLY in $cfg"; return 1; }
  edge=$((232 + (dly - 1100) / 2))
  export FADC_BASE_HI=$((edge - 32)) FADC_INT_LO=$edge FADC_INT_HI=$((edge + 184))
  WINDOWS_NOTE="FADC DLY $dly ns ($cfg) -> edge $edge: window_normal [0,$FADC_BASE_HI), integration [$FADC_INT_LO,$FADC_INT_HI)"
}

# IADC.  The pulse peaks at P = round((DLY - 1350) / 32.9) (runs 271-343,
# within 1 sample); the baseline is [0, P - 12), the pulse leaking in from
# P - 6, and the integration [P - 10, P + 22), >= 99.8 % of the mean pulse.
# P + 22 stays in a 126-sample record (RL 16) up to DLY ~4770.  Every board of
# the IADCT block must have the same DLY.
iadc_windows() {
  local cfg=$1 dlys dly peak
  [ -f "$cfg" ] || { echo "no run config $cfg for the IADC DLY"; return 1; }
  dlys=$(awk '/^[^ -]/ { in_b = ($1 == "IADCT:"); next } in_b' "$cfg" |
         grep -oE '\bDLY: *[0-9]+' | grep -oE '[0-9]+' | sort -u)
  [ "$(grep -c . <<<"$dlys")" -eq 1 ] || { echo "IADC DLY in $cfg is not one value: $(echo $dlys)"; return 1; }
  read -r dly peak IADC_BASE_HI IADC_INT_LO IADC_INT_HI < <(awk -v d="$dlys" 'BEGIN {
    d -= (d % 1000) % 16                  # the firmware drops the low 4 bits of the ns part
    p = int((d - 1350) / 32.9 + 0.5)
    print d, p, p - 12, p - 10, p + 22 }')
  export IADC_BASE_HI IADC_INT_LO IADC_INT_HI
  WINDOWS_NOTE="IADC DLY $dly ns ($cfg) -> peak $peak: window_normal [0,$IADC_BASE_HI), integration [$IADC_INT_LO,$IADC_INT_HI)"
}
