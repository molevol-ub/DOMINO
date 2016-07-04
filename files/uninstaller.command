#!/bin/sh
## Script for uninstalling and deleting files for DOMINO
current_dir=`pwd`
sh uninstall
rm -rf bin
rm -rf docs
rm -rf src
cd ..
rm -rf $current_dir
