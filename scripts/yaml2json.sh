#!/bin/bash
for f in `find . -name '*.yaml'`; do
    outdir=`dirname $f`
    python -c 'import sys, yaml, json; json.dump(yaml.load(sys.stdin), sys.stdout, indent=4)' < $f > ${outdir}/session_info.json
done

