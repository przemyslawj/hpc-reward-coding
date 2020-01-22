#!/bin/bash
# Usage ./rclone_filter_download <filter_name> <experiment_month>
set -x
rclone copy -P --filter-from $1 workdrive:cheeseboard/$2/ $2/
