#!/bin/bash
mkdir -p $1
rclone copy -P --filter-from filter_movie.txt prez_drive:cheeseboard/$1/ $1/
