#!/bin/bash

(return 0 2>/dev/null) && sourced=1 || sourced=0
if [ $sourced -eq 1 ]; then
   echo "

ci/clean.sh : This script is not destined to be sourced
             Please Type \`bash ci/clean.sh\` or \`ci/clean.sh\`

"
else

basedir=$(dirname $(realpath $0))
cd $basedir/..
$basedir/install.sh clean

fi
