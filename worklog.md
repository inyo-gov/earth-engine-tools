# Work Log

Use this log to capture decisions, data issues, and completion notes during
the WY24–25 update.

## Template

**Date:**
**Task:**
**Tags:**
**Inputs/Files:**
**Notes/Decisions:**
**Outputs:**
**Follow-ups:**

---

## Entries

**Date:** 2026-01-26  
**Task:** WY24–25 Earth Engine update setup and runs  
**Tags:** `ee`, `auth`, `zonal_stats`, `wy24-25`, `gitignore`  
**Inputs/Files:** `environment.yml`, `inyo_zonal_stats_col_2.ini`, `ee_zonal_stats_by_zone.py`, `README.md`, `.gitignore`  
**Notes/Decisions:** Created/updated conda env, upgraded Earth Engine API to support localhost auth, registered project and set `ee-icwd-2023`, switched INI to relative paths and WY24–25 dates, updated pandas compatibility in zonal stats, and archived legacy folders via `.gitignore`.  
**Outputs:** Runs completed for main 380 and c315; d75 completed; a_dss_200 started. Outputs in `stats*` (ignored).  
**Follow-ups:** Resume `a_dss_200` run if interrupted; move remaining legacy folders into `archive/`.  
