# ContactAngleCalculator cac
# dropletMolNumDensity unit number/A^3
# length unit: A
#command example：cac all.psf all.dcd xz 900 1000 20 {index 58158 to 68936} 3 0.0324 0.1 3 0.5 0.5
#command example：cac 10.xyz  10.xyz  xz 0 0 1 {index 60024 to 87023}  1 9.0  0.1  3  1.5  1.2

set tcl_precision 7

proc cac {file_1 file_2 type firstFrame lastFrame step dropletAtomRange atomNumberInMol \
          dropletMolNumDensity glRatio sphereFactor thicknessFactor sliceFactor} {

if {[string equal $file_1 $file_2]} {
	mol new $file_1 first $firstFrame last $lastFrame step $step waitfor all
} else {
	mol new $file_1
	mol addfile $file_2 first $firstFrame last $lastFrame step $step waitfor all
}

set logfile [open log.dat a+]
set caFile [open contact_angle.dat a+]
puts $logfile "  "
puts $logfile "  "
puts $logfile "********************** file: $file_1 $file_2 **********************"
puts $logfile "Command: cac $file_1 $file_2 $type $firstFrame $lastFrame $step {$dropletAtomRange}\
          $atomNumberInMol $dropletMolNumDensity $glRatio $sphereFactor $thicknessFactor $sliceFactor"
puts $caFile "  "
puts $caFile "********************** file: $file_1 $file_2 **********************"
puts $caFile "Command: cac $file_1 $file_2 $type $firstFrame $lastFrame $step {$dropletAtomRange}\
          $atomNumberInMol $dropletMolNumDensity $glRatio $sphereFactor $thicknessFactor $sliceFactor"
puts $caFile "frameID  solidPlane  layerRadius  h  contactAngle"


# preliminary variables, lists setting and parameters calculating----------------------------		  
set pi 3.14159
set molDistance [expr {pow(1.0/$dropletMolNumDensity,1.0/3.0)}]
set sphereFactor [expr {$sphereFactor*$molDistance}]
set thicknessFactor [expr {$thicknessFactor*$molDistance}]
set sliceFactor [expr {$sliceFactor*$molDistance}]
set threshold [expr {$glRatio*$dropletMolNumDensity*4.0*$pi*pow($sphereFactor,3)/3.0}]

set dropletAtomSet [atomselect top "$dropletAtomRange"]
set dropletAtomNum [$dropletAtomSet num]
set dropletAtomIndex [$dropletAtomSet get index]
set dropletMolIndex {}
for {set i 0} {$i < $dropletAtomNum} {incr i $atomNumberInMol} {
	lappend dropletMolIndex [lrange $dropletAtomIndex $i [expr {$i+$atomNumberInMol-1}]]	
}
set dropletMelNum [llength $dropletMolIndex]	

puts "   "
puts "    molDistance = $molDistance"
puts "         sphere = $sphereFactor"
puts "      threshold = $threshold"
puts "      thickness = $thicknessFactor"
puts "          slice = $sliceFactor"
puts $logfile "    molDistance = $molDistance"
puts $logfile "         sphere = $sphereFactor"
puts $logfile "      threshold = $threshold"
puts $logfile "      thickness = $thicknessFactor"
puts $logfile "          slice = $sliceFactor"
close $logfile
close $caFile

set n [molinfo top get numframes]
for {set i 0} {$i < $n} {incr i} {
	set logfile [open log.dat a+]
	set caFile [open contact_angle.dat a+]
	set frameID [expr {$i*$step+$firstFrame}]
	set dropletMolPosList {}
	set dropletAtomDelIndex {}
	set dropletAtomRemainIndex {}
	set dropletMolRemainPosList_X {}
	set dropletMolRemainPosList_Y {}
	set dropletMolRemainPosList_Z {}
	
	puts " "
	puts "---------------- frameID: $frameID ----------------"
	puts $logfile " "
	puts $logfile "---------------- frameID: $frameID ----------------"

# step_1: coarse-grain multi-atom molecules and remove the evaparated molecules-----------------------------
	foreach molAtomIndex $dropletMolIndex {	
		set molSel [atomselect top "index $molAtomIndex" frame $i]
		set molPos [measure center $molSel]	
		lappend dropletMolPosList $molPos
	}
	
	# loop every droplet molecules to detect if they are vapor molecules
	set molIndexOutCount 0
	set molDelNum 0
	set molRemainNum 0
	foreach molIndexOut $dropletMolPosList {
		set xOut [lindex $molIndexOut 0]
		set yOut [lindex $molIndexOut 1]
		set zOut [lindex $molIndexOut 2]
		
		set localLiquidMolNum 0
		foreach molIndexIn $dropletMolPosList {
			set xIn [lindex $molIndexIn 0]
			set yIn [lindex $molIndexIn 1]
			set zIn [lindex $molIndexIn 2]			
			set molDis [expr {pow(pow($xOut-$xIn,2.0)+pow($yOut-$yIn,2.0)+pow($zOut-$zIn,2.0),1.0/2.0)}]		
			if { $molDis < $sphereFactor } {
				set localLiquidMolNum [expr {$localLiquidMolNum+1}]
			}			
		}
		
		if { $localLiquidMolNum < $threshold} {	
			foreach molSplit [lindex $dropletMolIndex $molIndexOutCount] {
				lappend dropletAtomDelIndex $molSplit
			}
			set molDelNum [expr {$molDelNum+1}]
			set delMolAtomID [lindex $dropletMolIndex $molIndexOutCount]
			puts "Molecule Deleted: [expr {$molIndexOutCount+1}]/$dropletMelNum, atom index $delMolAtomID"
			puts $logfile "Molecule Deleted: [expr {$molIndexOutCount+1}]/$dropletMelNum, atom index $delMolAtomID"
		} else {
			foreach molSplit [lindex $dropletMolIndex $molIndexOutCount] {
				lappend dropletAtomRemainIndex $molSplit
			}
			lappend dropletMolRemainPosList_X [lindex $molIndexOut 0]
			lappend dropletMolRemainPosList_Y [lindex $molIndexOut 1]
			lappend dropletMolRemainPosList_Z [lindex $molIndexOut 2]
			set molRemainNum [expr {$molRemainNum+1}]
		}
		
		set molIndexOutCount [expr {$molIndexOutCount+1}]
	}

	set totalDropletMolNum [expr {$molRemainNum+$molDelNum}]
	if {$totalDropletMolNum != $dropletMelNum} {
		puts ""
		puts "****************************"
		puts "----------- frameID: $frameID -----------"
		puts "totalDropletMolNum != dropletMelNum!"
		puts "-------STOP HERE----------!"
		
		puts $logfile "totalDropletMolNum != dropletMelNum!"
		puts $logfile "-------STOP HERE----------!"
		puts "****************************"
		puts ""
		close $logfile
		~~~~~~~~~~HaHaHa~~~~~~~~~~
	}
	
	# output a PDB file of the remaining part	
	set delFlag [llength $dropletAtomDelIndex]
	if {$delFlag > 0} {
		set sysRemain [atomselect top "not index $dropletAtomDelIndex" frame $i]
		$sysRemain writepdb frame_$frameID.pdb
	} else {
			set sysRemain [atomselect top all frame $i]
			$sysRemain writepdb frame_$frameID.pdb
	}

# step_2: calculate the real contact area A, different wetting planes are considered
#         wetting on XZ plane
	if {[string equal -nocase $type xz]} {
		set dropletAtomSel [atomselect top "index $dropletAtomRemainIndex" frame $i]
		set dropletAtomRemainSize [measure minmax $dropletAtomSel]
		set xl [lindex $dropletAtomRemainSize 0 0]
		set xh [lindex $dropletAtomRemainSize 1 0]
		set zl [lindex $dropletAtomRemainSize 0 2]
		set zh [lindex $dropletAtomRemainSize 1 2]
		set substrateArea [atomselect top "not $dropletAtomRange and x>$xl and x<$xh and z>$zl and z<$zh" frame $i]
		set substrateSize [measure minmax $substrateArea]
		set yl [expr {[lindex $substrateSize 1 1]}]
		set yh [expr {$yl+$thicknessFactor}]
		
		
		set dropletMolBelow 0
		foreach molRemainIndex $dropletMolRemainPosList_Y {
			if {$molRemainIndex < $yl} {
				set dropletMolBelow [expr {$dropletMolBelow+1}]
			}
		}

		set dropletMolNumInCA [expr {$molRemainNum-$dropletMolBelow}]
		if {$delFlag > 0} {
			set dropletAtomInLayer [atomselect top "index $dropletAtomRemainIndex and y>$yl and y<$yh \
												and not index $dropletAtomDelIndex" frame $i]
		} else {
			set dropletAtomInLayer [atomselect top "index $dropletAtomRemainIndex and y>$yl and y<$yh" frame $i]
		}
												
		set dropletAtomInLayerIndex [$dropletAtomInLayer get index]
		set caFlag [llength $dropletAtomInLayerIndex]
		if {$caFlag > 1} {
			set dropletAtomInLayerSize [measure minmax $dropletAtomInLayer]
			set xl [expr {[lindex $dropletAtomInLayerSize 0 0]-0.001}]
			set xh [expr {[lindex $dropletAtomInLayerSize 1 0]+0.001}]
			set zl [expr {[lindex $dropletAtomInLayerSize 0 2]-0.001}]
			set zh [expr {[lindex $dropletAtomInLayerSize 1 2]+0.001}]
			set xSliceNumber [expr {int(($xh-$xl)/$sliceFactor)}]
			set xExtra [expr {($xh-$xl)-$sliceFactor*$xSliceNumber}]
			set zSliceNumber [expr int(($zh-$zl)/$sliceFactor)]
			set zExtra [expr {($zh-$zl)-$sliceFactor*$zSliceNumber}]

			puts " Wetting layer yl = $yl"
			puts "  dropletMolBelow = $dropletMolBelow"
			puts "dropletMolNumInCA = $dropletMolNumInCA"
			puts $logfile " Wetting layer yl = $yl"
			puts $logfile "  dropletMolBelow = $dropletMolBelow"
			puts $logfile "dropletMolNumInCA = $dropletMolNumInCA"

			set xContactArea 0
			for {set j 0} {$j<=$xSliceNumber} {incr j 1} {
				if {$j < $xSliceNumber} {
					set binWidth $sliceFactor
					set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$binWidth*$j}] [expr {$xl+$binWidth*($j+1)}] $yl $yh $zl $zh]
				} else {
					set binWidth $xExtra
					set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$sliceFactor*$j}] [expr {$xl+$sliceFactor*$j+$binWidth}] $yl $yh $zl $zh]
				}					
				set dropletAtomInBin [atomselect top "$bin and index $dropletAtomInLayerIndex" frame $i]
				set dropletAtomInBinSize [measure minmax $dropletAtomInBin]
				set zlBin [lindex $dropletAtomInBinSize 0 2]
				set zhBin [lindex $dropletAtomInBinSize 1 2]
				set xContactArea [expr {$xContactArea+$binWidth*($zhBin-$zlBin)}]
			}			
			set zContactArea 0
			for {set j 0} {$j<=$zSliceNumber} {incr j 1} {
				if {$j < $zSliceNumber} {
					set binWidth $sliceFactor
					set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						$xl $xh $yl $yh [expr {$zl+$binWidth*$j}] [expr {$zl+$binWidth*($j+1)}]]
				} else {
					set binWidth $zExtra
					set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						$xl $xh $yl $yh [expr {$zl+$sliceFactor*$j}] [expr {$zl+$sliceFactor*$j+$binWidth}]]
				}							
				set dropletAtomInBin [atomselect top "$bin and index $dropletAtomInLayerIndex" frame $i]
				set dropletAtomInBinSize [measure minmax $dropletAtomInBin]
				set xlBin [lindex $dropletAtomInBinSize 0 0]
				set xhBin [lindex $dropletAtomInBinSize 1 0]
				set zContactArea [expr {$zContactArea+$binWidth*($xhBin-$xlBin)}]
			}
				
	# step_3: construct an equivalent sessile droplet with volume V
			set xRadius [expr {sqrt($xContactArea/$pi)}]
			set zRadius [expr {sqrt($zContactArea/$pi)}]
			if {$xRadius < $zRadius} {
				set layerRadius $xRadius
			} else {
				set layerRadius $zRadius
			}
			puts "        xRadius = $xRadius"
			puts "        zRadius = $zRadius"
			puts "    layerRadius = $layerRadius"
			puts $logfile "          xRadius = $xRadius"
			puts $logfile "          zRadius = $zRadius"
			puts $logfile "      layerRadius = $layerRadius"
			
	# step_4: calculate the contact angle
			set vReal [expr {double($dropletMolNumInCA)/$dropletMolNumDensity}]
			set p [expr {3.0*$layerRadius*$layerRadius}]
			set q [expr {-6*$vReal/$pi}]
			set root [expr {sqrt($q*$q/4.0+$p*$p*$p/27.0)}]
			set h [expr {pow(-$q/2.0+$root,1.0/3.0)-pow($q/2.0+$root,1.0/3.0)}]	
			set CA [expr {2*180*atan(double($h)/$layerRadius)/$pi}]
			#set CA [expr {2*180*atan(2/1)/$pi}]
			puts "              h = $h"
			puts "             CA = $CA"
			puts $logfile "              h = $h"
			puts $logfile "             CA = $CA"
			puts $caFile [format "%d %.2f %.2f %.2f %.2f" $frameID $yl $layerRadius $h $CA]
		} else {
			puts "             CA = 180, no contact"
			puts $logfile "             CA = 180, No Contact"
			puts $caFile "$frameID 0.00 0.00 0.00 180.00"		
		}
		close $logfile
		close $caFile	
	
	} elseif {[string equal -nocase $type yz]} {
	
		set dropletAtomSel [atomselect top "index $dropletAtomRemainIndex" frame $i]
		set dropletAtomRemainSize [measure minmax $dropletAtomSel]
		set yl [lindex $dropletAtomRemainSize 0 1]
		set yh [lindex $dropletAtomRemainSize 1 1]
		set zl [lindex $dropletAtomRemainSize 0 2]
		set zh [lindex $dropletAtomRemainSize 1 2]
		set substrateArea [atomselect top "not $dropletAtomRange and y>$yl and y<$yh and z>$zl and z<$zh" frame $i]
		set substrateSize [measure minmax $substrateArea]
		set xl [expr {[lindex $substrateSize 1 0]}]
		set xh [expr {$xl+$thicknessFactor}]
		
		
		set dropletMolBelow 0
		foreach molRemainIndex $dropletMolRemainPosList_X {
			if {$molRemainIndex < $xl} {
				set dropletMolBelow [expr {$dropletMolBelow+1}]
			}
		}

		set dropletMolNumInCA [expr {$molRemainNum-$dropletMolBelow}]
		if {$delFlag > 0} {
			set dropletAtomInLayer [atomselect top "index $dropletAtomRemainIndex and x>$xl and x<$xh \
												and not index $dropletAtomDelIndex" frame $i]
		} else {
			set dropletAtomInLayer [atomselect top "index $dropletAtomRemainIndex and x>$xl and x<$xh" frame $i]
		}
												
		set dropletAtomInLayerIndex [$dropletAtomInLayer get index]
		set caFlag [llength $dropletAtomInLayerIndex]
		if {$caFlag > 1} {
			set dropletAtomInLayerSize [measure minmax $dropletAtomInLayer]
			set yl [expr {[lindex $dropletAtomInLayerSize 0 1]-0.001}]
			set yh [expr {[lindex $dropletAtomInLayerSize 1 1]+0.001}]
			set zl [expr {[lindex $dropletAtomInLayerSize 0 2]-0.001}]
			set zh [expr {[lindex $dropletAtomInLayerSize 1 2]+0.001}]
			set ySliceNumber [expr {int(($yh-$yl)/$sliceFactor)}]
			set yExtra [expr {($yh-$yl)-$sliceFactor*$ySliceNumber}]
			set zSliceNumber [expr int(($zh-$zl)/$sliceFactor)]
			set zExtra [expr {($zh-$zl)-$sliceFactor*$zSliceNumber}]

			puts " Wetting layer xl = $xl"
			puts "  dropletMolBelow = $dropletMolBelow"
			puts "dropletMolNumInCA = $dropletMolNumInCA"
			puts $logfile " Wetting layer xl = $xl"
			puts $logfile "  dropletMolBelow = $dropletMolBelow"
			puts $logfile "dropletMolNumInCA = $dropletMolNumInCA"

			set yContactArea 0
			for {set j 0} {$j<=$ySliceNumber} {incr j 1} {
				if {$j < $ySliceNumber} {
					set binWidth $sliceFactor
					set bin [format "y>%.3f and y<%.3f and x>%.3f and x<%.3f and z>%.3f and z<%.3f" \
						[expr {$yl+$binWidth*$j}] [expr {$yl+$binWidth*($j+1)}] $xl $xh $zl $zh]
				} else {
					set binWidth $yExtra
					set bin [format "y>%.3f and y<%.3f and x>%.3f and x<%.3f and z>%.3f and z<%.3f" \
						[expr {$yl+$sliceFactor*$j}] [expr {$yl+$sliceFactor*$j+$binWidth}] $xl $xh $zl $zh]
				}					
				set dropletAtomInBin [atomselect top "$bin and index $dropletAtomInLayerIndex" frame $i]
				set dropletAtomInBinSize [measure minmax $dropletAtomInBin]
				set zlBin [lindex $dropletAtomInBinSize 0 2]
				set zhBin [lindex $dropletAtomInBinSize 1 2]
				set yContactArea [expr {$yContactArea+$binWidth*($zhBin-$zlBin)}]
			}			
			set zContactArea 0
			for {set j 0} {$j<=$zSliceNumber} {incr j 1} {
				if {$j < $zSliceNumber} {
					set binWidth $sliceFactor
					set bin [format "y>%.3f and y<%.3f and x>%.3f and x<%.3f and z>%.3f and z<%.3f" \
						$yl $yh $xl $xh [expr {$zl+$binWidth*$j}] [expr {$zl+$binWidth*($j+1)}]]
				} else {
					set binWidth $zExtra
					set bin [format "y>%.3f and y<%.3f and x>%.3f and x<%.3f and z>%.3f and z<%.3f" \
						$yl $yh $xl $xh [expr {$zl+$sliceFactor*$j}] [expr {$zl+$sliceFactor*$j+$binWidth}]]
				}							
				set dropletAtomInBin [atomselect top "$bin and index $dropletAtomInLayerIndex" frame $i]
				set dropletAtomInBinSize [measure minmax $dropletAtomInBin]
				set ylBin [lindex $dropletAtomInBinSize 0 1]
				set yhBin [lindex $dropletAtomInBinSize 1 1]
				set zContactArea [expr {$zContactArea+$binWidth*($yhBin-$ylBin)}]
			}
				
	# step_3: construct an equivalent sessile droplet with volume V
			set yRadius [expr {sqrt($yContactArea/$pi)}]
			set zRadius [expr {sqrt($zContactArea/$pi)}]
			if {$yRadius < $zRadius} {
				set layerRadius $yRadius
			} else {
				set layerRadius $zRadius
			}
			puts "        yRadius = $yRadius"
			puts "        zRadius = $zRadius"
			puts "    layerRadius = $layerRadius"
			puts $logfile "          yRadius = $yRadius"
			puts $logfile "          zRadius = $zRadius"
			puts $logfile "      layerRadius = $layerRadius"
			
	# step_4: calculate the contact angle
			set vReal [expr {double($dropletMolNumInCA)/$dropletMolNumDensity}]
			set p [expr {3.0*$layerRadius*$layerRadius}]
			set q [expr {-6*$vReal/$pi}]
			set root [expr {sqrt($q*$q/4.0+$p*$p*$p/27.0)}]
			set h [expr {pow(-$q/2.0+$root,1.0/3.0)-pow($q/2.0+$root,1.0/3.0)}]	
			set CA [expr {2*180*atan(double($h)/$layerRadius)/$pi}]
			#set CA [expr {2*180*atan(2/1)/$pi}]
			puts "              h = $h"
			puts "             CA = $CA"
			puts $logfile "              h = $h"
			puts $logfile "             CA = $CA"
			puts $caFile [format "%d %.2f %.2f %.2f %.2f" $frameID $xl $layerRadius $h $CA]
		} else {
			puts "             CA = 180, no contact"
			puts $logfile "             CA = 180, No Contact"
			puts $caFile "$frameID 0.00 0.00 0.00 180.00"		
		}
		close $logfile
		close $caFile	
	
	} elseif {[string equal -nocase $type xy]} {
	
		set dropletAtomSel [atomselect top "index $dropletAtomRemainIndex" frame $i]
		set dropletAtomRemainSize [measure minmax $dropletAtomSel]
		set xl [lindex $dropletAtomRemainSize 0 0]
		set xh [lindex $dropletAtomRemainSize 1 0]
		set yl [lindex $dropletAtomRemainSize 0 1]
		set yh [lindex $dropletAtomRemainSize 1 1]
		set substrateArea [atomselect top "not $dropletAtomRange and x>$xl and x<$xh and y>$yl and y<$yh" frame $i]
		set substrateSize [measure minmax $substrateArea]
		set zl [expr {[lindex $substrateSize 1 2]}]
		set zh [expr {$zl+$thicknessFactor}]
		
		
		set dropletMolBelow 0
		foreach molRemainIndex $dropletMolRemainPosList_Z {
			if {$molRemainIndex < $zl} {
				set dropletMolBelow [expr {$dropletMolBelow+1}]
			}
		}

		set dropletMolNumInCA [expr {$molRemainNum-$dropletMolBelow}]
		if {$delFlag > 0} {
			set dropletAtomInLayer [atomselect top "index $dropletAtomRemainIndex and z>$zl and z<$zh \
												and not index $dropletAtomDelIndex" frame $i]
		} else {
			set dropletAtomInLayer [atomselect top "index $dropletAtomRemainIndex and z>$zl and z<$zh" frame $i]
		}
												
		set dropletAtomInLayerIndex [$dropletAtomInLayer get index]
		set caFlag [llength $dropletAtomInLayerIndex]
		if {$caFlag > 1} {
			set dropletAtomInLayerSize [measure minmax $dropletAtomInLayer]
			set xl [expr {[lindex $dropletAtomInLayerSize 0 0]-0.001}]
			set xh [expr {[lindex $dropletAtomInLayerSize 1 0]+0.001}]
			set yl [expr {[lindex $dropletAtomInLayerSize 0 1]-0.001}]
			set yh [expr {[lindex $dropletAtomInLayerSize 1 1]+0.001}]
			set xSliceNumber [expr {int(($xh-$xl)/$sliceFactor)}]
			set xExtra [expr {($xh-$xl)-$sliceFactor*$xSliceNumber}]
			set ySliceNumber [expr int(($yh-$yl)/$sliceFactor)]
			set yExtra [expr {($yh-$yl)-$sliceFactor*$ySliceNumber}]

			puts " Wetting layer zl = $zl"
			puts "  dropletMolBelow = $dropletMolBelow"
			puts "dropletMolNumInCA = $dropletMolNumInCA"
			puts $logfile " Wetting layer zl = $zl"
			puts $logfile "  dropletMolBelow = $dropletMolBelow"
			puts $logfile "dropletMolNumInCA = $dropletMolNumInCA"

			set xContactArea 0
			for {set j 0} {$j<=$xSliceNumber} {incr j 1} {
				if {$j < $xSliceNumber} {
					set binWidth $sliceFactor
					set bin [format "x>%.3f and x<%.3f and z>%.3f and z<%.3f and y>%.3f and y<%.3f" \
						[expr {$xl+$binWidth*$j}] [expr {$xl+$binWidth*($j+1)}] $zl $zh $yl $yh]
				} else {
					set binWidth $xExtra
					set bin [format "x>%.3f and x<%.3f and z>%.3f and z<%.3f and y>%.3f and y<%.3f" \
						[expr {$xl+$sliceFactor*$j}] [expr {$xl+$sliceFactor*$j+$binWidth}] $zl $zh $yl $yh]
				}					
				set dropletAtomInBin [atomselect top "$bin and index $dropletAtomInLayerIndex" frame $i]
				set dropletAtomInBinSize [measure minmax $dropletAtomInBin]
				set ylBin [lindex $dropletAtomInBinSize 0 1]
				set yhBin [lindex $dropletAtomInBinSize 1 1]
				set xContactArea [expr {$xContactArea+$binWidth*($yhBin-$ylBin)}]
			}			
			set yContactArea 0
			for {set j 0} {$j<=$ySliceNumber} {incr j 1} {
				if {$j < $ySliceNumber} {
					set binWidth $sliceFactor
					set bin [format "x>%.3f and x<%.3f and z>%.3f and z<%.3f and y>%.3f and y<%.3f" \
						$xl $xh $zl $zh [expr {$yl+$binWidth*$j}] [expr {$yl+$binWidth*($j+1)}]]
				} else {
					set binWidth $yExtra
					set bin [format "x>%.3f and x<%.3f and z>%.3f and z<%.3f and y>%.3f and y<%.3f" \
						$xl $xh $zl $zh [expr {$yl+$sliceFactor*$j}] [expr {$yl+$sliceFactor*$j+$binWidth}]]
				}							
				set dropletAtomInBin [atomselect top "$bin and index $dropletAtomInLayerIndex" frame $i]
				set dropletAtomInBinSize [measure minmax $dropletAtomInBin]
				set xlBin [lindex $dropletAtomInBinSize 0 0]
				set xhBin [lindex $dropletAtomInBinSize 1 0]
				set yContactArea [expr {$yContactArea+$binWidth*($xhBin-$xlBin)}]
			}
				
	# step_3: construct an equivalent sessile droplet with volume V
			set xRadius [expr {sqrt($xContactArea/$pi)}]
			set yRadius [expr {sqrt($yContactArea/$pi)}]
			if {$xRadius < $yRadius} {
				set layerRadius $xRadius
			} else {
				set layerRadius $yRadius
			}
			puts "        xRadius = $xRadius"
			puts "        yRadius = $yRadius"
			puts "    layerRadius = $layerRadius"
			puts $logfile "          xRadius = $xRadius"
			puts $logfile "          yRadius = $yRadius"
			puts $logfile "      layerRadius = $layerRadius"
			
	# step_4: calculate the contact angle
			set vReal [expr {double($dropletMolNumInCA)/$dropletMolNumDensity}]
			set p [expr {3.0*$layerRadius*$layerRadius}]
			set q [expr {-6*$vReal/$pi}]
			set root [expr {sqrt($q*$q/4.0+$p*$p*$p/27.0)}]
			set h [expr {pow(-$q/2.0+$root,1.0/3.0)-pow($q/2.0+$root,1.0/3.0)}]	
			set CA [expr {2*180*atan(double($h)/$layerRadius)/$pi}]
			#set CA [expr {2*180*atan(2/1)/$pi}]
			puts "              h = $h"
			puts "             CA = $CA"
			puts $logfile "              h = $h"
			puts $logfile "             CA = $CA"
			puts $caFile [format "%d %.2f %.2f %.2f %.2f" $frameID $zl $layerRadius $h $CA]
		} else {
			puts "             CA = 180, no contact"
			puts $logfile "             CA = 180, No Contact"
			puts $caFile "$frameID 0.00 0.00 0.00 180.00"		
		}
		close $logfile
		close $caFile	
	}
}	
}
