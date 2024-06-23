#!/bin/csh -f

## Distribute jobs to multiple cores

# Arguments
if ($#argv < 3) then
    echo "$0 n_core i_core flconfig.out"
    exit 2
endif
set ncore=$argv[1]
set icore=$argv[2]
set fldata=$argv[3]

# User variables
set template=./flsim_conf
set outhead=./flsim_out
set atmhead=./flsim_atm
set imghead=./img_flsim
set dathead=./les05_
set a3dtail=.mdla3d
set phstail=.mdlphs
set sfcfile=" "
set nph=5e8
set ws=(0045 0055 0067)
set ss=(0) # RT solver options: F3D and ICA
@ njob = $#ws * $#ss

# System variables
set tmp0=tmp.$$.0
set mcarats=../bin/mcarats

# Get configuration parameters for this experiment
set buf = (`tail -n +2 $fldata | wc -l`)
set nt = $buf[1]
set xpos = (`tail -n +2 $fldata | awk '{print $2}'`)
set ypos = (`tail -n +2 $fldata | awk '{print $3}'`)
set zloc = (`tail -n +2 $fldata | awk '{print $4}'`)
set the1 = (`tail -n +2 $fldata | awk '{print $5}'`)
set phi1 = (`tail -n +2 $fldata | awk '{print $6}'`)
set the0 = (`tail -n +2 $fldata | awk '{print $7}'`)
set phi0 = (`tail -n +2 $fldata | awk '{print $8}'`)
set psi1 = (`tail -n +2 $fldata | awk '{print $9}'`)

# Loop for all jobs
@ ijob = $icore
while($ijob <= $njob)
    @ iw = ( $ijob - 1 ) / $#ss + 1
    @ is = $ijob - ( $iw - 1 ) * $#ss
    set s=$ss[$is]
    set w=$ws[$iw]
    echo "#### (iw, is) = "$iw $is

    # Make a namelist file
    rm -f $tmp0
    grep -A 10000 mcarWld_nml_init $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Wld_njob = $nt >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarSca_nml_init $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Sca_inpfile = \'$dathead$w$phstail\' >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarAtm_nml_init $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Atm_inpfile = \'$dathead$w$a3dtail\' >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarSfc_nml_init $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Sfc_inpfile = \'$sfcfile\' >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarSrc_nml_init $template | grep -B 10000 -m 1 ^/ >> $tmp0
    grep -A 10000 mcarFlx_nml_init $template | grep -B 10000 -m 1 ^/ >> $tmp0
    grep -A 10000 mcarRad_nml_init $template | grep -B 10000 -m 1 ^/ >> $tmp0
    grep -A 10000 mcarVis_nml_init $template | grep -B 10000 -m 1 ^/ >> $tmp0
    grep -A 10000 mcarPho_nml_init $template | grep -B 10000 -m 1 ^/ >> $tmp0
    grep -A 10000 mcarWld_nml_job  $template | grep -B 10000 -m 1 ^/ >> $tmp0
    grep -A 10000 mcarAtm_nml_job  $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    cat $atmhead"_wav"$w >> $tmp0
    echo Atm_idread = 1 >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarSfc_nml_job  $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Sfc_idread = 1 >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarSrc_nml_job  $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Src_the = `echo $the0[1] | awk '{print 180-$1}'` >> $tmp0
    echo Src_phi = `echo $phi0[1] | awk '{print 180+$1}'` >> $tmp0
    echo / >> $tmp0
    grep -A 10000 mcarRad_nml_job  $template | grep -B 10000 -m 1 ^/ | grep -v / >> $tmp0
    echo Rad_xpos = $xpos[1] >> $tmp0
    echo Rad_ypos = $ypos[1] >> $tmp0
    echo Rad_zloc = $zloc[1] >> $tmp0
    echo Rad_the  = $the1[1] >> $tmp0
    echo Rad_phi  = $phi1[1] >> $tmp0
    echo Rad_psi  = $psi1[1] >> $tmp0
    echo / >> $tmp0
    @ it = 2
    while($it <= $nt)
	echo "!----------------" >> $tmp0
	echo \&mcarWld_nml_job / >> $tmp0
	echo \&mcarAtm_nml_job Atm_idread = 0 / >> $tmp0
	echo \&mcarSfc_nml_job Sfc_idread = 0 / >> $tmp0
	echo \&mcarSrc_nml_job >> $tmp0
	echo Src_the = `echo $the0[$it] | awk '{print 180-$1}'` >> $tmp0
	echo Src_phi = `echo $phi0[$it] | awk '{print 180+$1}'` >> $tmp0
	echo / >> $tmp0
	echo \&mcarRad_nml_job >> $tmp0
	echo Rad_xpos = $xpos[$it] >> $tmp0
	echo Rad_ypos = $ypos[$it] >> $tmp0
	echo Rad_zloc = $zloc[$it] >> $tmp0
	echo Rad_the  = $the1[$it] >> $tmp0
	echo Rad_phi  = $phi1[$it] >> $tmp0
	echo Rad_psi  = $psi1[$it] >> $tmp0
	echo / >> $tmp0
	@ it = $it + 1
    end
    #cat $tmp0

    # Run mcarats
    set outfile=$outhead"_s"$s"_wav"$w
    echo Running mcarats for $outfile"..."
    $mcarats $nph $s $tmp0 $outfile

    @ ijob = $ijob + $ncore
end

rm -f $tmp0
