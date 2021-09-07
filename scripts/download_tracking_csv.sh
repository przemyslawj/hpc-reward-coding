#!/bin/bash
rclone copy -P --filter-from filter_tracking_csv.txt workdrive:cheeseboard/$1/ $1/
