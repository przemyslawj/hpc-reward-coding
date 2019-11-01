#!/bin/bash
# Usage ./rclone_filter_download <filter_name> <experiment_month>
set -x
rclone copy -P --filter-from $1 workdrive:cheeseboard-down/down_2/$2/ $2/
