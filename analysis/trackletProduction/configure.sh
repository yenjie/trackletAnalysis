if [ $# -lt 1 ]; then
    echo "Usage: ./configure.sh [5/8] [0/1/2]"
    echo "       0: DATA, 1: EPOS, 2: HIJING"
    exit 1
fi

if [ $1 -eq 5 -a $2 -eq 2 ]; then
    echo "  ! 5TeV HIJING sample not available"
    exit 1
fi

case $1 in
    5)
        _PILEUP=0.003;;
    8)
        _PILEUP=0.006;;
    *)
        echo "5/8 only"
        exit 1;;
esac

case $2 in
    0)
        _SAMPLE=DATA;;
    1)
        _SAMPLE=EPOS;;
    2)
        _SAMPLE=HIJING;;
    *)
        echo "0/1/2 only"
        exit 1;;
esac

echo "resetting #defines..."
sed -i '/#define _CONFIG_START/,/#define _CONFIG_END/{//!d}' definitions.h

echo -e "setting defaults for ${_SAMPLE} at ${1}TeV...\n"
case $2 in
    0)
        cat definitions.h;;
    [12])
        awk '/#define _CONFIG_START/ {
            print;
            print "#define _PILEUP '${_PILEUP}'";
            print "#define _DATA_VTX_'${1}'TEV";
            print "#define _'${_SAMPLE}'_VTX_'${1}'TEV";
            print "#define _'${_SAMPLE}'_MULTWGHT_'${1}'TEV";
            print "#define _POKE_HOLES";
            print "#define _HOLES_'${1}'TEV"; next;
        }1' definitions.h | tee definitions.tmp && mv definitions.tmp definitions.h
        echo -e "\n"
        ;;
    *)
        echo "!!!"
esac
