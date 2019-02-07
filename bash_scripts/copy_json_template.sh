#!/bin/bash
#
## ucdbioinfo_supernova_pipeline copy_json_template.sh
## small script to create a template project setup json file to stdout

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

call="cat $DIR/../templates/run_setup.json"
eval $call
