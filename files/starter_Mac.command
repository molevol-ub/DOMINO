#!/bin/sh
## starter
#import variables
current_dir=$( cd $(dirname $0) ; pwd -P )
export QT_QPA_FONTDIR=$current_dir"/"fonts 
$current_dir"/"DOMINO.app"/"Contents"/"MacOS"/"DOMINO
