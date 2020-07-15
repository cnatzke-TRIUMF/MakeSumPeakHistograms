#!/bin/bash

ANALDIR=/tig/belodon_data4/S9038/S9038-Nov2019/AnalysisTrees
HISTDIR=/tig/belodon_data4/S9038/Projects/SumPeakAnalysis/MakeHistograms/myOutput
SCALDIR=/tig/belodon_data4/S9038/Projects/SumPeakAlignment/myOutput/LinearCalibrationParameters/UnformattedOutput

DLEN=${#ANALDIR}
SORTCODE=/tig/belodon_data4/S9038/Projects/SumPeakAnalysis/MakeHistograms/myBuild/SumPeakHistograms

# calibration files
CALFILE=/tig/belodon_data4/S9038/S9038-Nov2019/PreSrcCal.cal
#CALFILE=/tig/belodon_data4/S9038/S9038-Nov2019/PostSrcCal.cal
#CFILE=/tig/belodon_data4/S9038/PostSrcCal.cal

m=15600 #Only Sort after this run First good run: 52042
n=99999 #Only Sort until this run
BAD_RUN_LIST=(15719 15720)

#############################################
# Helper functions
#############################################

# checks if item is in list
bad_run_check() {
   [[ $1 =~ (^|[[:space:]])$2($|[[:space:]]) ]] && return 0 || return 1
}

#############################################
# Main
#############################################

if [ $# -eq 0 ]
then
   echo "Pass 1 or 2 arguments"
   echo "1: $0 run_number"
   echo "   Process all subruns in run_number"
   echo "2: $0 run_number batch_size"
   echo "   Process all subruns in run_number in batches of batch_size"

   #for f in $ANALDIR/*.root
#
   #do
   #   #echo "$f"
   #   g=${f:DLEN+9}
   #   h=${g:0:${#g}-6}
   #   i=${g:0:${#g}-9} # run number
#
   #   RUNFILES=$ANALDIR/analysis${i}_"*".root
   #   HFILE=$HISTDIR/TwoPhotonGriffinAnglePlots${i}_"*".root
#
   #   if [ $i -gt $m ] && [ $n -gt $i ];
   #   then
   #      if [ ! -f $HFILE ];
   #      then
   #         echo "$SORTCODE $RUNFILES $SELECTOR --max-workers=4"
   #         #$SORTCODE $RUNFILES $SELECTOR --max-workers=4
   #         mv ./TwoPhotonGriffinAnglePlots* $HISTDIR/
   #      fi
   #      if [ -f $HFILE ]
   #      then
   #         # if the sort code gets updated or the analysis tree is updated
   #         if [ $HFILE -ot $SORTCODE ] || [ $HFILE -ot $f ];
   #         then
   #            echo "$SORTCODE $RUNFILES $SELECTOR --max-workers=4"
   #            #$SORTCODE $RUNFILES $SELECTOR --max-workers=4
   #            mv TwoPhotonGriffinAnglePlots* $HISTDIR/
   #         fi
   #      fi
   #   fi
   #done
fi

# for a single run
if [ $# -eq 1 ]
then
    file=$ANALDIR/analysis"$1"_*.root;
    g=${file:DLEN+9}
    h=${g:0:${#g}-6}
    i=${g:0:${#g}-7} # run number

      RUNFILES=$ANALDIR/analysis${i}_"*".root
      HFILE=$HISTDIR/gg_${i}.root

      if [ $i -gt $m ] && [ $n -gt $i ];
      then
         if [ ! -f $HFILE ];
         then
           # echo "$SORTCODE $RUNFILES $CALFILE"
           # echo "mv gg_* $HISTDIR/"
           $SORTCODE $RUNFILES $CALFILE
           mv gg_* $HISTDIR/gg_${i}_total.root

         fi
         if [ -f $HFILE ]
         then
            # if the sort code gets updated or the analysis tree is updated
            if [ $HFILE -ot $SORTCODE ] || [ $HFILE -ot $f ];
            then
               echo "Sort code or analysis tree have changed, rerunning sort"
                echo "$SORTCODE $RUNFILES $CALFILE"
                echo "mv gg_* $HISTDIR/"
               #$SORTCODE $RUNFILES $SELECTOR --max-workers=4
               #mv TwoPhotonGriffinAnglePlots* $HISTDIR/
            fi
         fi
     elif [[ $i < $m ]]; then
         echo "Run $1 is less than the low run cutoff: $m"
         echo "Exiting ..."
     else
         echo "$1 is non accepted run format"
         echo "Exiting ..."
      fi
fi

# for a range of subruns
if [ $# -eq 2 ] # Run Number, Batch size
then
    declare -a FILE_BATCH
    declare -a FILE_LIST
    FILE_ITER=1
    BATCH_COUNTER=0

    RUNFILES=$ANALDIR/analysis${i}_"*".root
    HFILE=$HISTDIR/energy_angle_${i}.root
    SCALFILE=$SCALDIR/secondary_linear_cal_parameters_${1}.txt

    if [ ! -f "$SCALFILE" ]; then
        echo "::: ERROR ::: Cannot find"
        echo "$SCALFILE"
        echo "::: EXITING ::: "
        exit 1
    fi

    FILE_LIST=$(find ${ANALDIR}/ -name analysis"$1"_???.root | sort)
    NUM_FOUND=$(find ${ANALDIR}/ -name analysis"$1"_???.root | wc -l)
    TOTAL_BATCHES=$((NUM_FOUND / $2 + 1))
    echo "Found ${NUM_FOUND} Files"
    echo "Requires ${TOTAL_BATCHES} Batches"

    for file in $FILE_LIST;
    do

        FILE_BATCH=("${FILE_BATCH[@]}" "$file")

        # Process files in specified batch sizes
        if (( $FILE_ITER % $2 == 0 ))
        then
            #echo "Files: ${FILE_BATCH[@]}"

            if [ $1 -gt $m ] && [ $n -gt $1 ];
            then
                if [ ! -f $HFILE ];
                then
                    #echo "$SORTCODE ${FILE_BATCH[@]} $CALFILE $SCALFILE"
                    #echo "mv secondary_calibrated_histograms_* $HISTDIR/energy_angle_${1}_${BATCH_COUNTER}.root"
                    $SORTCODE ${FILE_BATCH[@]} $CALFILE $SCALFILE
                    mv secondary_calibrated_histograms* $HISTDIR/energy_angle_${1}_${BATCH_COUNTER}.root
                fi
                if [ -f $HFILE ]
                then
                   # if the sort code gets updated or the analysis tree is updated
                    if [ $HFILE -ot $SORTCODE ] || [ $HFILE -ot $file ];
                    then
                        echo "Sort code or analysis tree have changed, rerunning sort"
                        echo "$SORTCODE $FILE_BATCH $CALFILE"
                        echo "mv gg_* $HISTDIR/gg_$1_$BATCH_COUNTER.root"
                    fi
                fi
            fi

            unset FILE_BATCH
            BATCH_COUNTER=$((BATCH_COUNTER + 1))
        fi

        # Process remaining files
        if (( $BATCH_COUNTER == $((TOTAL_BATCHES - 1))  && (( $FILE_ITER == $NUM_FOUND)) ))
        then
            #echo "Files: ${FILE_BATCH[@]}"
            if [ $1 -gt $m ] && [ $n -gt $1 ];
            then
                if [ ! -f $HFILE ];
                then
                    #echo "$SORTCODE ${FILE_BATCH[@]} $CALFILE $SCALFILE"
                    #echo "mv secondary_calibrated_histograms_* $HISTDIR/energy_angle_${1}_${BATCH_COUNTER}.root"
                    $SORTCODE ${FILE_BATCH[@]} $CALFILE $SCALFILE
                    mv secondary_calibrated_histograms* $HISTDIR/energy_angle_${1}_${BATCH_COUNTER}.root
                fi
                if [ -f $HFILE ]
                then
                   # if the sort code gets updated or the analysis tree is updated
                    if [ $HFILE -ot $SORTCODE ] || [ $HFILE -ot $file ];
                    then
                        echo "Sort code or analysis tree have changed, rerunning sort"
                        echo "$SORTCODE $FILE_BATCH $CALFILE"
                        echo "mv gg_* $HISTDIR/gg_$1_$BATCH_COUNTER.root"
                    fi
                fi
            fi

            unset FILE_BATCH
        fi

        FILE_ITER=$((FILE_ITER + 1))

    done
fi

if [[ $# > 2 ]]
then
   echo "Pass 1 or 2 arguments"
   echo "1: $0 run_number"
   echo "   Process all subruns in run_number"
   echo "2: $0 run_number batch_size"
   echo "   Process all subruns in run_number in batches of batch_size"
fi
