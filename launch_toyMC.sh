N_EV=$1
N_J=$2

N_EV_PER_J=$(($N_EV/$N_J))

rm cmds
for i in $(seq 1 $N_J); do
echo "./toyMC_main \"_pt_small\" ${N_EV_PER_J} ${i}" >> cmds
done

g++ -O3 toyMC_main.cxx -o toyMC_main `root-config --cflags --glibs`

cat cmds | xargs -P ${N_J} -I CMD bash -c CMD
