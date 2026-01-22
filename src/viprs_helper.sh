#!/usr/bin/env

make_viprs_sumstats() {

    # ---- arguments ----
    local piece1=$1
    local piece2=$2
    local piece3=$3
    local path_outputs=$4

    # ---- basic checks ----
    if [[ $# -ne 4 ]]; then
        echo "Usage: make_viprs_sumstats <piece1> <piece2> <piece3> <output_dir>" >&2
        return 1
    fi

    if [[ ! -d "${path_outputs}" ]]; then
        mkdir -p "${path_outputs}" || return 1
    fi

    # ---- piece 1 ----
    awk '
    BEGIN {
        OFS="\t";
        print "#CHROM","POS","ID","REF","ALT1","A1","A1_FREQ","OBS_CT","BETA","SE","T_STAT","P"
    }
    NR>1 && $6>0 && $7>0 && $7<=1 {
        chrom=$1
        pos=$3
        id=$2
        a1=toupper($4)
        beta=$5
        se=$6
        p=$7
        tstat=beta/se

        print chrom,pos,id,"NA",a1,a1,"NA","NA",beta,se,tstat,p
    }' "${piece1}" > "${path_outputs}/piece1.viprs.tsv"

    # ---- piece 2 ----
    awk '
    BEGIN {
        OFS="\t";
        print "#CHROM","POS","ID","REF","ALT1","A1","A1_FREQ","OBS_CT","BETA","SE","T_STAT","P"
    }
    NR>1 && $6>0 {
        id=$1
        chrom=$2
        a1=toupper($3)
        ref=toupper($4)
        beta=$5
        se=$6
        n=$7
        tstat=beta/se

        print chrom,"NA",id,ref,a1,a1,"NA",n,beta,se,tstat,"NA"
    }' "${piece2}" > "${path_outputs}/piece2.viprs.tsv"

    # ---- piece 3 ----
    awk '
    BEGIN {
        OFS="\t";
        print "#CHROM","POS","ID","REF","ALT1","A1","A1_FREQ","OBS_CT","BETA","SE","T_STAT","P"
    }
    NR>1 && $5=="ADD" && $7!=0 && $9>0 && $9<=1 {
        chrom=$1
        pos=$3
        id=$2
        a1=toupper($4)
        beta=$7
        tstat=$8
        p=$9
        n=$6
        se=beta/tstat

        if (se>0)
            print chrom,pos,id,"NA",a1,a1,"NA",n,beta,se,tstat,p
    }' "${piece3}" > "${path_outputs}/piece3.viprs.tsv"

    # ---- split by chromosome ----
    pushd "${path_outputs}" >/dev/null || return 1

    for f in piece*.viprs.tsv; do
        awk '
        NR==1 {header=$0; next}
        {
            out = FILENAME ".chr" $1
            print >> out
        }
        END {
            for (i in seen)
                print header > i
        }' "$f"
    done

    popd >/dev/null
}
