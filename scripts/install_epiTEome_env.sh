#### install the damn epiTEome env
conda create -y -n epiTEome_env -c conda-forge mamba
conda activate epiTEome_env

mamba install -y -c conda-forge python=3.11 perl=5.32
mamba install -y -c bioconda bedtools=2.30 segemehl=0.2.0 perl-app-cpanminus perl-bioperl-core perl-bio-samtools perl-set-intervaltree perl-statistics-descriptive

# pip install ngsutils==0.5.7
# # change all py2 files into py3 <<<as the epiTEome developer said: "Die Ficker!">>>
# 2to3 -w /home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils/

###############################################################################
#                                                                             #
# use samtools <<<updated - version 1.16>>> instead of the 'ngsutils' package #
#                                                                             #
# samtools >= 1.10 is needed!                                                 #
# check if there is '-N' argument in the filtering section                    #
#                                                                             #
###############################################################################

#------------------------------------------

### main directory
mkdir /home/yoyerush/yo/methylome_pipeline/transposition_epiTEome/
cd /home/yoyerush/yo/methylome_pipeline/transposition_epiTEome/
wget https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/mto1_paper_R_scripts/epiTEome_te_insertion_sites_script.sh
chmod +x ./epiTEome_te_insertion_sites_script.sh

### scripts
mkdir epiteome_scripts
cd epiteome_scripts
wget https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/epiTEome_Yo_edit.pl # edited for using 'samtool=1.16' instead of the 'ngsutils' package
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/epiTEome.pl # 'epiTEome.pl' original script
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/idxEpiTEome.pl
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/insertTEsintoFasta.pl
chmod +x ./epiTEome_Yo_edit.pl
chmod +x ./epiTEome.pl
chmod +x ./idxEpiTEome.pl
chmod +x ./insertTEsintoFasta.pl
cd ../

### genome files
mkdir genome_indx
cd genome_indx
wget -O TAIR10_chr_all.fna.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
gunzip TAIR10_chr_all.fna.gz
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/test/tair10TEs.gff3
cd ../

### test data
mkdir test_data
cd test_data
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/test/teid.lst
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/test/Chr2.fasta
wget https://raw.githubusercontent.com/jdaron/epiTEome/master/test/unmapped.fastq.bz2
cd ../

#------------------------------------------

### create target list (TE IDs list)
mkdir te_lists
cd te_lists
# run to create retro-TEs list:
awk '$9 ~ /sF=LTR\/Copia|sF=LTR\/Gypsy|sF=LINE\/L1/ {
     split($9, arr, ";");
     for (i in arr) {
         if (arr[i] ~ /ID=/) {
             split(arr[i], id, "=");
             print id[2];
         }
     }
}' ../genome_indx/tair10TEs.gff3 | sort | uniq > retrotransposons_list.txt

# run to create helitron - "ATREP9" family list:
awk '$9 ~ /fam=ATREP9/ {
     split($9, arr, ";");
     for (i in arr) {
         if (arr[i] ~ /ID=/) {
             split(arr[i], id, "=");
             print id[2];
         }
     }
}' ../genome_indx/tair10TEs.gff3 | sort | uniq > helitron_atrep9_list.txt

# run to create copia superfamily list:
awk '$9 ~ /sF=LTR\/Copia/ {
     split($9, arr, ";");
     for (i in arr) {
         if (arr[i] ~ /ID=/) {
             split(arr[i], id, "=");
             print id[2];
         }
     }
}' ../genome_indx/tair10TEs.gff3 | sort | uniq > copia_list.txt

# run to create gypsy superfamily list:
awk '$9 ~ /sF=LTR\/Gypsy/ {
     split($9, arr, ";");
     for (i in arr) {
         if (arr[i] ~ /ID=/) {
             split(arr[i], id, "=");
             print id[2];
         }
     }
}' ../genome_indx/tair10TEs.gff3 | sort | uniq > gypsy_list.txt

# run to create LINE/L1 superfamily list:
awk '$9 ~ /sF=LINE\/L1/ {
     split($9, arr, ";");
     for (i in arr) {
         if (arr[i] ~ /ID=/) {
             split(arr[i], id, "=");
             print id[2];
         }
     }
}' ../genome_indx/tair10TEs.gff3 | sort | uniq > line_l1_list.txt

cd ../

#------------------------------------------

###########  another option:
###########  don't use this anymore! but if using 'ngsutils', this could help...

###################################
## if there is error writing: << Unknown command 'filter' >>
#  /home/yoyerush/anaconda3/envs/epiTEome_env/bin/bamutils
#  or find it with the command << which bamutils >>
#  note: also to fastqutils
#
#  then make this modifications inside the bamutil file:

#  #!/bin/sh
#  REAL=/home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils
#  DIR=`dirname "$REAL"`/../
#  SUBDIR=/home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils/bam/
# 
#  usage() {
#      echo "Usage: $(basename $0) COMMAND"
#      echo ""
#      cat $DIR/ngsutils/$SUBDIR/README
#      echo ""
#      echo "Run '$(basename $0) help CMD' for more information about a specific command"
#      echo -n "ngsutils "
# 
#      cd $DIR
#      GV="$(git show master --format='%h %ai' | head -n 1)"
#      VERSION=$(echo "$GV" | awk '{print $1}')
#      echo "$(cat VERSION | sed -e 's/\n//')-$VERSION"
# 
#      exit 1
#  }
# 
#  if [ "$1" = "" ]; then
#      usage
#  fi
# 
# 
#  . "/home/yoyerush/anaconda3/pkgs/python-3.7.6-h0371630_2/lib/python3.7/venv/scripts/common/activate"
# 
#  export PYTHONPATH=$PYTHONPATH:"$DIR"
# 
#  if [ "$SUBDIR" = "ngs" ]; then
#      if [ -e "$DIR"/.git ]; then
#          if [ "$1" = "update" ]; then
#              cd "$DIR"
# 
#              echo "Updating from current branch"
#              git pull
# 
#              exit 0
#          elif [ "$1" = "switch" ]; then
#              cd "$DIR"
# 
#              if [ "$2" != "" ]; then
#                  echo "Switching to branch: $2"
#                  git fetch origin
#                  if [ "$(git branch -l | grep "$2")" = "" ]; then
#                      git branch -t $2 origin/$2
#                  fi
#                  git checkout $2
#                  git pull
#              else
#                  git branch -a
#              fi
#              exit 0
#          fi
#      fi
#  fi
# 
# 
#  if [ "$1" = "help" ]; then
#      if [ "$2" = "" ]; then
#          usage
#      fi
# 
#      action=/home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils/bam/filter.py
# 
#      if [ ! -e "$action" ]; then
#          action="$DIR"/ngsutils/$SUBDIR/$2.sh
#          if [ ! -e "$action" ]; then
#              echo "Unknown command '$2'"
#              exit 1
#          fi
#      fi
#      "$action" -h
# 
#  elif [ "$1" = "version" ]; then
#      cd $DIR
#      GV="$(git show master --format='%h %ai' | head -n 1)"
#      VERSION=$(echo "$GV" | awk '{print $1}')
#      echo "$(cat VERSION | sed -e 's/\n//')-$VERSION"
# 
#  elif [ "$1" = "profile" ]; then
#      shift
#      action=/home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils/bam/filter.py
# 
#      if [ ! -e "$action" ]; then
#          action="$DIR"/ngsutils/$SUBDIR/filter.sh
#          if [ ! -e "$action" ]; then
#              echo "Unknown command 'filter'"
#              exit 1
#          fi
#      fi
#      shift
# 
#  #    i=0
#  #    for arg in "$@"; do
#  #        ARGS[$i]="$arg"
#  #        ((++i))
#  #    done
# 
#      echo "Saving profile information to profile.output" 1>&2
#      exec python -m cProfile -o profile.output "$DIR"/ngsutils/$SUBDIR/$action "$@"
#  else
#      action=/home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils/bam/filter.py
# 
#      if [ ! -e "$action" ]; then
#          action="$DIR"/ngsutils/$SUBDIR/filter.sh
#          if [ ! -e "$action" ]; then
#              echo "Unknown command '$1'"
#              exit 1
#          fi
#      fi
#      shift
# 
#      exec "$action" "$@"
#  fi


#### edit fastqutils
#  #!/bin/bash
#  REAL=/home/yoyerush/anaconda3/envs/epiTEome_env/lib/python3.11/site-packages/ngsutils
#  DIR=`dirname "$REAL"`/../
#  SUBDIR=$(basename $0 | sed -e 's/utils//')
#  
#  function usage() {
#      echo "Usage: $(basename $0) COMMAND"
#      echo ""
#      cat $DIR/ngsutils/$SUBDIR/README
#      echo ""
#      echo "Run '$(basename $0) help CMD' for more information about a specific command"
#      echo -n "ngsutils "
#      
#      cd $DIR
#      GV="$(git show master --format='%h %ai' | head -n 1)"
#      VERSION=$(echo "$GV" | awk '{print $1}')
#      echo "$(cat VERSION | sed -e 's/\n//')-$VERSION"
#  
#      exit 1
#  }
#  
#  if [[ "$1" == "" ]]; then
#      usage
#  fi
#  
#  
#  . "/home/yoyerush/anaconda3/pkgs/python-3.7.6-h0371630_2/lib/python3.7/venv/scripts/common/activate"
#  
#  export PYTHONPATH=$PYTHONPATH:"$DIR"
#  
#  if [[ "$SUBDIR" == "ngs" ]]; then
#      if [[ -e "$DIR"/.git ]]; then
#          if [[ "$1" == "update" ]]; then
#              cd "$DIR"
#      
#              echo "Updating from current branch"
#              git pull 
#  
#              exit 0
#          elif [[ "$1" == "switch" ]]; then
#              cd "$DIR"
#      
#              if [[ "$2" != "" ]]; then
#                  echo "Switching to branch: $2"
#                  git fetch origin
#                  if [[ "$(git branch -l | grep "$2")" == "" ]]; then
#                      git branch -t $2 origin/$2
#                  fi
#                  git checkout $2
#                  git pull
#              else
#                  git branch -a
#              fi
#              exit 0
#          fi
#      fi
#  fi
#  
#  
#  if [[ "$1" == "help" ]]; then
#      if [[ "$2" == "" ]]; then
#          usage
#      fi
#      
#      action=$2.py
#      
#      if [[ ! -e "$DIR"/ngsutils/$SUBDIR/$action ]]; then
#          action=$2.sh
#          if [[ ! -e "$DIR"/ngsutils/$SUBDIR/$action ]]; then
#              echo "Unknown command '$2'"
#              exit 1
#          fi
#      fi
#      "$DIR"/ngsutils/$SUBDIR/$action -h
#  
#  elif [[ "$1" == "version" ]]; then
#      cd $DIR
#      GV="$(git show master --format='%h %ai' | head -n 1)"
#      VERSION=$(echo "$GV" | awk '{print $1}')
#      echo "$(cat VERSION | sed -e 's/\n//')-$VERSION"
#  
#  elif [[ "$1" == "profile" ]]; then
#      shift
#      action=$1.py
#  
#      if [[ ! -e "$DIR"/ngsutils/$SUBDIR/$action ]]; then
#          action=$1.sh
#          if [[ ! -e "$DIR"/ngsutils/$SUBDIR/$action ]]; then
#              echo "Unknown command '$1'"
#              exit 1
#          fi
#      fi
#      shift
#  
#      ARGS=()
#      i=0
#      for arg in "$@"; do
#          ARGS[$i]="$arg"
#          ((++i))
#      done
#      
#      echo "Saving profile information to profile.output" 1>&2
#      exec python -m cProfile -o profile.output "$DIR"/ngsutils/$SUBDIR/$action "${ARGS[@]}"
#  else
#      action=$1.py
#  
#      if [[ ! -e "$DIR"/ngsutils/$SUBDIR/$action ]]; then
#          action=$1.sh
#          if [[ ! -e "$DIR"/ngsutils/$SUBDIR/$action ]]; then
#              echo "Unknown command '$1'"
#              exit 1
#          fi
#      fi
#      shift
#  
#      ARGS=()
#      i=0
#      for arg in "$@"; do
#          ARGS[$i]="$arg"
#          ((++i))
#      done
#      
#      exec "$DIR"/ngsutils/$SUBDIR/$action "${ARGS[@]}"
#  fi
