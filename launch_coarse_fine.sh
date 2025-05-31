#compile program
g++ -O3 fit_coarse_fine.cxx -o fit `root-config --glibs --cflags`

#create commands
rm cmds
for i in $(seq 0 9); do
  j=$(($i+1))
  echo "./fit ${i} ${j} \"foo_dat_0\"" >> cmds
done

#launch commands
cat cmds | xargs -P 10 -I CMD bash -c CMD
