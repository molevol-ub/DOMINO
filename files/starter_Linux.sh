#!/bin/sh
## starter
#import variables
current_dir=$( cd $(dirname $0) ; pwd -P )

echo $current_dir
export QT_QPA_FONTDIR=$current_dir"/"fonts 
export LD_LIBRARY_PATH=$current_dir"/"lib

$current_dir"/"DOMINO
