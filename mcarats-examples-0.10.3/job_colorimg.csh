#!/bin/csh -f

## To generate color images from binary (GrADS compatible)  data files

# Arguments
if ($#argv < 6) then
    echo "$0 blue.ctl green.ctl red.ctl outfileHead pwr rmax (timeLag timeWid fmax)"
    exit 2
endif
set imgc=./$argv[4]
set pwr=$argv[5]
set rmax=$argv[6]

# User variables
set fcol=(1.28 1.0 0.8) # color balance

# System variables
set bin_exposure=../bin/bin_exposure
set bin_gray=../bin/bin_gray
set tmp0=./$$.tmp.0
set imgs=(./$$.tmp_img_B ./$$.tmp_img_G ./$$.tmp_img_R)
set imgt=./$$.tmp_img.tif


# Gray images
if ($rmax == 0) then # automatic exposure
    echo $bin_exposure $argv[7] $argv[8] $argv[9] $argv[5] $argv[1] $argv[2] $argv[3]
    $bin_exposure $argv[7] $argv[8] $argv[9] $argv[5] $argv[1] $argv[2] $argv[3] > $tmp0
    foreach c(1 2 3)
	$bin_gray $fcol[$c] $rmax $pwr $argv[$c] $imgs[$c] < $tmp0
    end
else # manual exposure
    foreach c(1 2 3)
	$bin_gray $fcol[$c] $rmax $pwr $argv[$c] $imgs[$c]
    end
endif
set nx=(`grep XDEF $argv[1] | awk '{print $2}'`)
set ny=(`grep YDEF $argv[1] | awk '{print $2}'`)
set nz=(`grep ZDEF $argv[1] | awk '{print $2}'`)
set nt=(`grep TDEF $argv[1] | awk '{print $2}'`)
set nv=(`grep VARS $argv[1] | awk '{print $2}'`)

# Color images
@ t = 1
while($t <= $nt)
    if ($t <= 9) then
	set tstr="000"$t
    else if ($t <= 99) then
	set tstr="00"$t
    else if ($t <= 999) then
	set tstr="0"$t
    else
	set tstr=$t
    endif
    @ v = 1
    while($v <= $nv)
	@ z = 1
	while($z <= $nz)
	    if ( -r $imgs[1]"_t"$t"_var"$v"_z"$z".gray" ) then
		set str="_t"$t"_var"$v"_z"$z
		set str1="_t"$tstr"_var"$v"_z"$z
		set pwr1=$pwr
		@ ntry = 0
		while($ntry < 30) # loop for trials
		    #// This is because the gray image file is sometimes broken.
		    set err=0
		    foreach c(1 2 3)
			convert -depth 8 -size $nx"x"$ny+0 $imgs[$c]$str".gray" $imgs[$c]$str".png"
			if ($status != 0) then
			    set err=1
			    echo Error: $imgs[$c]$str".gray" "with pwr1="$pwr1
			endif
		    end
		    if ($err == 0) then # success
			set ntry=1000000
		    else # error, then retry
			set pwr1=`echo $pwr1"*1.00333" | bc`
			foreach c(1 2 3)
			    if ($rmax == 0) then # automatic exposure
				$bin_gray $fcol[$c] $rmax $pwr1 $argv[$c] $imgs[$c] $z $v $t < $tmp0
			    else # manual exposure
				$bin_gray $fcol[$c] $rmax $pwr1 $argv[$c] $imgs[$c] $z $v $t
			    endif
			end
		    endif
		    @ ntry = $ntry + 1
		end
		composite -compose CopyGreen $imgs[2]$str".png" $imgs[3]$str".png" $imgt
		composite -compose CopyBlue  $imgs[1]$str".png" $imgt $imgc$str1".tif"
		echo $imgc$str1".tif" : complete
		rm -f $imgs[1]$str* $imgs[2]$str* $imgs[3]$str* $imgt
	    endif
	    @ z = $z + 1
	end
	@ v = $v + 1
    end
    @ t = $t + 1
end
rm -f $tmp0
