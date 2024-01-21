#!/bin/bash
woh_file=$1
wh_file=$2
log_file=$3

if (( ${#woh_file} <= 0)); then echo "Error: Please provide an input pdb file."; exit 0; fi
if [ ! -f $woh_file ]; then echo "Error: Input file \"$1\" does not exist."; exit 0; fi
if (( ${#wh_file} <= 0)); then echo "Error: Please provide an output filename."; exit 0; fi
if (( ${#log_file} <= 0)); then log_file="/dev/null"; fi
if ! command -v tleap &> /dev/null; then echo "Error: The Amber program \"tleap\" is required for this script."; exit 0; fi
if ! command -v python &> /dev/null; then echo "Error: Python is required for this script."; exit 0; fi

woh_file_temp="$wh_file.temp"
while [ -f $woh_file_temp ]; do woh_file_temp="$woh_file_temp.temp"; done

sed '/HOH/d' $woh_file > $woh_file_temp

tleap -f > $log_file - <<_EOF
source leaprc.protein.ff14SB
mol = loadpdb $woh_file_temp
savepdb mol $wh_file
quit
_EOF

#rm leap.log

sed -i 's/HIE/HIS/g' $wh_file

read -r -d '' python_command << EOM
import fileinput
data = []
last_id = ""
for line in fileinput.input('$woh_file_temp'):
    if line[:4] == 'ATOM' and last_id != line[22:26]:
        last_id = line[22:26]
        data.append(line[21])
count = -1
last_id = ""
for line in fileinput.input('$wh_file', inplace=True):
    if line[:4] == 'ATOM' and last_id != line[22:26]:
        last_id = line[22:26]
        count += 1
    if line[:4] == 'ATOM':
        try:
            print(line[:21] + data[count] + line[22:], end='')
        except:
            print(line, end='')  
    else:
        print(line, end='')
EOM

python -c "$python_command"

rm "$woh_file_temp"
