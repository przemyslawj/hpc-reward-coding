#!/bin/bash
# Usage ./upload_down2 <filter_name> <experiment_month>
set -x
rclone copy -P --filter-from $1 $2/ workdrive:cheeseboard-down/down_2/$2/
