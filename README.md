# earth-engine-tools
Earth Engine Python tools for extracting Landsat and GRIDMET zonal statistics
over groundwater dependent ecosystem parcels.

## Setup
- Create the conda env: `conda env create -f environment.yml`
- Activate it: `conda activate ee-tools`
- Authenticate EE (local browser): `earthengine authenticate --auth_mode=localhost`
- Set project (required): `earthengine set_project <your-ee-project-id>`

If you get a "project is not registered" error, register the project in the
Earth Engine console before proceeding.

## Yearly Update Workflow
Use `inyo_zonal_stats_col_2.ini` as the control file (relative paths are
expected). Current zone inputs live in `gis/`; the old `shapefiles/` folder
is archived.

1) Update date range in `[INPUTS]`:
   - `start_year = 2024`
   - `end_year = 2025`
2) Choose the dataset by setting `zone_shp_path` and `output_workspace`:
   - Example: `zone_shp_path = gis/parcels_rasterized.shp`
   - Example: `output_workspace = stats`
3) Run the zonal stats script:
   - `python ee_zonal_stats_by_zone.py -i inyo_zonal_stats_col_2.ini`

Repeat step 2–3 for each dataset if there are multiple sets, different project, etc. Runs
can take hours; do them one at a time.

## Outputs
Outputs are written to the `output_workspace` directory and organized by
zone name in subfolders (e.g., `stats/<ZONE>/...csv`).

The script also writes derived files next to the input shapefile
(`.geojson` and `_tiles.json`) if they don't already exist.

These outputs are large and are intentionally ignored by git (see `.gitignore`).
