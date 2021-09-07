#!/bin/bash
rclone copy -P --filter-from filter_dlc_result.txt drive_synced:cheeseboard/$1/ $1/
