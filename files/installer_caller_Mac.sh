#!/bin/sh
## Open a terminal and show output of installing process
## $1 installdir
chmod 777 installer_Mac.command
script_sent="installer_Mac.command"
command="open -a Terminal.app installer_Mac.command"
exec $command & pid_sh=$!
while [[ $(ps | grep $script_sent) ]]
do
	sleep 30
done