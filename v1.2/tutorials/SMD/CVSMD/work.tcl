#!/bin/tclsh

##########################################################
### PULLING FORCES INTEGRATION TCL SCRIPT FOR TINKER-HP ##
##########################################################
#
## Owner: FREDERIC CELERSE
## Institution: SORBONNE UNIVERSITE
## Last update: 29/10/2019
#

#### Open the log file for reading and the output .dat file for writing
set file [open SMD_output1.dat r]
set output [open work.dat w]
##### Parameters input from user.
set nx -0.628		
set ny 0.717
set nz -0.302
set v 0.01
set dt 0.010
set outputfreq 100
set initdist 0
##### DO NOT MOVE FROM HERE !
set work 0.0
while { [gets $file line] != -1 } {
##### Determine if a line contains "SMD ". If so, write the
##### timestep followed by the calculated pulling work from 
##### the output file
        if {[lindex $line 0] == "SMD"} {
		set work [expr $work + $v*$dt*$outputfreq*($nx*[lindex $line 5] + $ny*[lindex $line 6] + $nz*[lindex $line 7])]
		puts $output "$initdist $work"	
		set initdist [expr $initdist + $v*$dt*$outputfreq]
        }
}
##### Close the log file and the output .dat file
close $file
close $output

                                




