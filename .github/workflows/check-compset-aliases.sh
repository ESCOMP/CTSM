#!/usr/bin/env bash
set -e

# Check that compset aliases return matching longnames

# Change to top level of clone
cd "$(git rev-parse --show-toplevel)"

# Check that query_config can run without error
cime/scripts/query_config --compsets 1>/dev/null

# Save previous IFS line-splitting behavior to restore at the end (keep this at the beining)
OLD_IFS=$IFS
IFS=$'\n'

bad_compsets() {
  # Subroutine to Find bad compsets
  # Set input arguments
  local ALIAS="$1"
  local LCOMPSET="$2"
  set +e
  # Relies on case sensitivity here: Alias should have $ALIAS and longname should have $LCOMPSET
  # --intervert-match can be shortened to -v
  # --fixed-strings can be shortened to -F
  local bad_compsets="$(cime/scripts/query_config --compsets | sort | uniq | grep --fixed-strings -- \"$ALIAS\" | grep --fixed-strings --invert-match -- \"$LCOMPSET\")"
  set -e
  if [[ "${bad_compsets}" != "" ]]; then
      echo "One or more compsets with $ALIAS alias but not $LCOMPSET longname:" >&2
      echo $bad_compsets  >&2
      exit 1
  fi
}

check_missing_alias_matches() {
  # Subroutine to Find bad compsets by checking that if something is not found in the alias -- then something else SHOULD be found in the long-compset name
  # Allows for regular expressions in the alias argument
  # Set input arguments
  local ALIAS="$1"
  local LCOMPSET="$2"
  set +e
  # Relies on case sensitivity here: Alias should have $ALIAS and longname should have $LCOMPSET
  # --intervert-match can be shortened to -v
  # --fixed-strings can be shortened to -F
  local bad_compsets="$(cime/scripts/query_config --compsets clm | awk 'NR > 5{print $0}' | grep -E --invert-match -- \"$ALIAS\" | grep --fixed-strings -- \"$LCOMPSET\")"
  set -e
  if [[ "${bad_compsets}" != "" ]]; then
      echo "One or more compsets without $ALIAS alias but not $LCOMPSET longname:" >&2
      echo $bad_compsets  >&2
      exit 1
  fi
}

# Now call the subroutine for various options
# -- Physics versions ---
bad_compsets Clm60 CLM60
bad_compsets Clm50 CLM50
bad_compsets Ctsm50 "CLM50%NWP"
# -- stub ROF or GLC
bad_compsets Rs SROF
bad_compsets Gs SGLC
# -- General veg type --
bad_compsets Bgc "%BGC"
bad_compsets BgcCrop "%BGC-CROP"
bad_compsets Clm50Sp "%SP"
bad_compsets Clm60Sp "%SP"
# ---- Note that specifying SP for NWP compsets as they shouldn't use BGC ----
bad_compsets Nwp "%NWP-SP"
# ---- NWP cases should always be stug GLC ---
bad_compsets Nwp Gs
bad_compsets Fates "%FATES"
bad_compsets FatesSp "%FATES-SP"
# ---- Using full NoAnthro This fails with the error: grep: thro ": No such file or directory
bad_compsets 'NoA' 'NOANTHRO_'

# -- Period of simulation --
bad_compsets I2000 2000_
bad_compsets I2010 2010_
bad_compsets I1850 1850_
bad_compsets IHist HIST_
bad_compsets ISSP585 SSP585_
bad_compsets ISSP126 SSP126_
bad_compsets ISSP119 SSP119_
bad_compsets ISSP245 SSP245_
bad_compsets ISSP370 SSP370_
bad_compsets ISSP460 SSP460_
bad_compsets ISSP534 SSP534_
# --- DATM forcing types ---
bad_compsets I1Pt "DATM%1PT"
bad_compsets Crujra "DATM%CRUJRA2024b"
bad_compsets Gswp "DATM%GSWP3v1"
bad_compsets "Cru " "DATM%CRU"
bad_compsets Nldas "DATM%NLDAS2"
bad_compsets Qian '_DATM%QIA_'
# ---- Using full "Spinup" This fails with the error: grep: up ": No such file or directory
bad_compsets 'Spi' 'DATM%CPLHIST'

# --- ROF models ---
bad_compsets Miz _MIZUROUTE_
bad_compsets "Rtm " _RTM_
bad_compsets RtmFl "_RTM%FLOOD_"

# --- CISM options ---
bad_compsets "G " '_CISM2%GRIS-EVOLVE_'
bad_compsets "Ga " '_CISM2%AIS-EVOLVE_'
bad_compsets "Gag " '_CISM2%AIS-EVOLVE%GRIS-EVOLVE_'


#
# Now check compsets that don't have something in the alias to make sure something else is there
#
check_missing_alias_matches "Rs|Rtm|Miz" "_MOSART_"
check_missing_alias_matches "Gs|G |Gag|Ga " "_DGLC%NOEVOLVE_"
check_missing_alias_matches 'Crujra|Cru|Nldas|As|1Pt|Qia|Spi' 'GSWP3v1'

# Restore line-splitting behavior (keep this at the end)
IFS=$OLD_IFS

exit 0