# Requirements:
* [R datatrace](https://github.com/przemyslawj/datatrace/) package installed

# Usage
## Merging tracking data and cell traces:
1. Prepare scripts for downloading the data.
Enter your data directory and copy there scripts for downloading the traces:
`cp cheeseboard_analysis/scripts/* .`
Set `EXP_MONTH` variable, for example:
`export EXP_MONTH=2019-08`

2. Download traces with filtered components if not there already:
`./download_down2.sh filter_msresults.txt $EXP_MONTH`

3. Download tracking csv files:
`./download_tracking_csv.sh $EXP_MONTH`

3. Prepare `session_info.yaml` files
`./yaml2json.sh`

4. Run matlab scripts merging the data:
From cheeseboard_analysis/matlab directory run matlab and start the run
`runCaimanForExp.m` script

5. Upload merged data traces csv to GDrive
`./upload_down2.sh filter_csv_traces.txt $EXP_MONTH`

## Data traces analysis scripts
- `active_cells.Rmd`: Change in the number of active cells and their activity
between days
- `create_place_field_df.R`: Calculates place fields and their statistics and writes them to an environment file.
- `place_field_stats.Rmd`: Analysis of place fields extracted with the script above.
- `goal_mutual_info.Rmd`: Calculates mutual information between variables marking the presence at the goal and the cell activity.
- `position_decoder.Rmd`: Bayes decoding of position given the cell activity.

The notebooks use metadata files about the trials which can be downloaded with
GDrive:
`./download.sh filter_trials_meta.txt $EXP_MONTH`
