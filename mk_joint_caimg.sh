#!/bin/bash

OUTPUT_DIR="$(pwd)/joined_caimg"
echo $OUTPUT_DIR
mkdir -p $OUTPUT_DIR

declare -A session_counter
for DATED_DIR in `ls -d 2019-*/`; do
  for SUBJECT_DIR in `ls -d $DATED_DIR/results_caimg/*`; do
      SUBJECT=$(basename $SUBJECT_DIR)
      SESSION_ID=${session_counter["${SUBJECT}"]}
      if [ ! ${SESSION_ID} ]; then
          SESSION_ID=1
      fi

      if [ ! -f ${SUBJECT_DIR} ]; then
          mkdir -p "${OUTPUT_DIR}/${SUBJECT}"
          JOINT_DIR=${OUTPUT_DIR}/${SUBJECT}/jointExtraction
          mkdir -p ${JOINT_DIR}
          mkdir -p ${JOINT_DIR}/alignment
          mkdir -p ${JOINT_DIR}/extracted
          mkdir -p ${JOINT_DIR}/sorted

      fi
      for SESSION_DIR in `ls -d ${SUBJECT_DIR}/Session*/`; do
          SESSION=$(basename $SESSION_DIR)
          DEST_DIR=${OUTPUT_DIR}/${SUBJECT}/Session${SESSION_ID}
          SRC_DIR=${SESSION_DIR}
          echo "Session ${SESSION_ID} is linked to: ${SRC_DIR}"
          ln -sf $(pwd)/${SRC_DIR} ${DEST_DIR}
          SESSION_ID=$((SESSION_ID + 1))
      done
      session_counter["${SUBJECT}"]=${SESSION_ID}
  done
done

ln -sf $(pwd)/session_info_template.csv ${OUTPUT_DIR}/session_info_template.csv

