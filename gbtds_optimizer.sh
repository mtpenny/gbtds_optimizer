#!/bin/bash -x

if [ $# -ne 2 ]; then
    echo "Usage: $0 <map_commands> <layout>"
    exit
fi

map_commands=$(cat $1)
map_file=$(echo $map_commands | awk 'BEGIN{mapopt=-1}{for(i=1;i<=NF;i++){if($(i)=="--map-filename" || $(i)=="-m") mapopt=i+1; if(i==mapopt) print $(i)}}')
maproot=${1%.*}
layout=$2

nftmp=${layout#layout_}
nfields=${nftmp%%f*}

echo map_file $map_file

#--map-filename fidu_mass6_rate.yield.csv --map-cadence 14.7315 --map-texp 42.56 --fields-filename

#Field dimensions
read lmin lcent lmax bmin bcent bmax < <(awk 'BEGIN{lmin=1e50; lmax=-1e50; bmin=1e50; bmax=-1e50}$4==0&&NR>1{if($2<lmin) lmin=$2; if($2>lmax) lmax=$2; if($3<bmin) bmin=$3; if($3>bmax) bmax=$3;}END{print lmin,0.5*(lmin+lmax),lmax,bmin,0.5*(bmin+bmax),bmax}' field_layouts/${nfields}fields/$layout.centers)

read padlm padlp padbm padbp < <(awk 'BEGIN{lmin=1e50; lmax=-1e50; bmin=1e50; bmax=-1e50}NF==3{if($2<lmin) lmin=$2; if($2>lmax) lmax=$2; if($3<bmin) bmin=$3; if($3>bmax) bmax=$3;}END{print lmin,lmax,bmin,bmax}' sca_layout.txt)

read maplmin maplmax mapbmin mapbmax < <(awk -v FS=',' 'BEGIN{lmin=1e50; lmax=-1e50; bmin=1e50; bmax=-1e50}NR>1{if($1<lmin) lmin=$1; if($1>lmax) lmax=$1; if($2<bmin) bmin=$2; if($2>bmax) bmax=$2;}END{print lmin,lmax,bmin,bmax}' $map_file)

step=0.2

lrange=$(echo 1 | awk -v lmin=$lmin -v lmax=$lmax -v lcent=$lcent -v padlm=$padlm -v padlp=$padlp -v maplmin=$maplmin -v maplmax=$maplmax -v step=$step 'function floor(x){if(x<0){if(int(x)==x) return int(x); else return int(x)-1}else return int(x)}{maxf=maplmax-padlp-(lmax-lcent); minf=maplmin-padlm-(lmin-lcent); print step*floor(maxf/step),-step*floor(-minf/step)}')

brange=$(echo 1 | awk -v bmin=$bmin -v bmax=$bmax -v bcent=$bcent -v padbm=$padbm -v padbp=$padbp -v mapbmin=$mapbmin -v mapbmax=$mapbmax -v step=$step 'function floor(x){if(x<0){if(int(x)==x) return int(x); else return int(x)-1}else return int(x)}{maxf=mapbmax-padbp-(bmax-bcent); minf=mapbmin-padbm-(bmin-bcent); print step*floor(maxf/step),-step*floor(-minf/step)}')

echo lmin lcent lmax bmin bcent bmax
echo $lmin $lcent $lmax $bmin $bcent $bmax
echo padlm padlp padbm padbp
echo $padlm $padlp $padbm $padbp
echo maplmin maplmax mapbmin mapbmax
echo $maplmin $maplmax $mapbmin $mapbmax
echo lrange
echo $lrange
echo brange
echo $brange


#lrange="5 -5"
#brange="-3 3"


python gbtds_optimizer.py $map_commands \
       --fields-filename field_layouts/${nfields}fields/${layout}.centers \
       --lrange $lrange --brange $brange \
       --lstep $step --bstep $step \
       --cadence-bounds 7.0 12.0 --nread-bounds 10 40 \
       --output-root ${maproot}.${layout}

python results_plotter.py -i $maproot.$layout --contour-resolution 5 \
       --smoothing 0.3 --lrange $maplmax $maplmin --brange $mapbmin $mapbmax \
       --save png --save-root $maproot.$layout --no-show

