for iapo in apo holo hK181A fV48A fQ123A fR5A
do
 for iadj in henm_hessian_adj
 do
   itmp=1.00
   ftmp=''  
    for isrc in  050 104 130 225
    do
      for isnk in 304 337 431 433
      do
	echo $iapo $iadj $isrc $isnk  

        sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_source.dat > input.dat
	time ./paths_TWISP.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_source_3D.dat > input.dat
	time ./paths_TWISP_3D.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill.dat > input.dat
	time ./paths_distill_multi.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_3D.dat > input.dat
	time ./paths_distill_3D_multi.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze.dat > input.dat
	time ./paths_analyze_multi.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_3D.dat > input.dat
	time ./paths_analyze_3D_multi.e

        sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_source_broken.dat > input.dat
	time ./paths_TWISP_broken.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_source_broken_3D.dat > input.dat
	time ./paths_TWISP_broken_3D.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_broken.dat > input.dat
	time ./paths_distill_broken_multi.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_broken_3D.dat > input.dat
	time ./paths_distill_broken_3D_multi.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_broken.dat > input.dat
	time ./paths_analyze_broken_multi.e
	sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/AAA/'$isrc'/g' -e 's/ZZZ/'$isnk'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_broken_3D.dat > input.dat
 	time ./paths_analyze_broken_3D_multi.e

      done
    done
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_all.dat > input.dat
    time ./paths_distill_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_all.dat > input.dat
    time ./paths_analyze_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_deg.dat > input.dat
    time ./paths_deg_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_partition.dat > input.dat
    time ./paths_partition.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_3D_all.dat > input.dat
    time ./paths_distill_3D_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_3D_all.dat > input.dat
    time ./paths_analyze_3D_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_deg_3D.dat > input.dat
    time ./paths_deg_3D_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_partition_3D.dat > input.dat
    time ./paths_partition_3D.e

    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_broken_all.dat > input.dat
    time ./paths_distill_broken_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_broken_all.dat > input.dat
    time ./paths_analyze_broken_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_deg_broken.dat > input.dat
    time ./paths_deg_broken_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_partition_broken.dat > input.dat
    time ./paths_partition_broken.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_distill_broken_3D_all.dat > input.dat
    time ./paths_distill_broken_3D_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_analyze_broken_3D_all.dat > input.dat
    time ./paths_analyze_broken_3D_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_deg_broken_3D.dat > input.dat
    time ./paths_deg_broken_3D_multi.e
    sed -e 's/YYY/'$iapo'/g' -e 's/XXX/'$iadj'/g' -e 's/TTT/'$itmp'/g' -e 's/WWW/'$ftmp'/g' input0_partition_broken_3D.dat > input.dat
    time ./paths_partition_broken_3D.e
#    
   time python edge_effect_multi.py $iapo ${iadj::`expr ${#iadj} - 4`}
   time python edge_effect_node_multi.py $iapo ${iadj::`expr ${#iadj} - 4`}
  done
done
