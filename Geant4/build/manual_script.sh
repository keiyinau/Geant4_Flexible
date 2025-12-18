./sim run_t20_100k_Co60.mac >> output_Co60.log
mv output_0.root Co60_Nearest.root
./sim run_t20_100k_Cs137.mac >> output_Cs137.log
mv output_0.root Cs137_Nearest.root
./sim run_t20_100k.mac >> output_Na22.log
mv output_0.root Na22_nearest.root
#./sim cosmic_muon.mac >> output_cosmic.log
#mv output_0.root Cosmic_nearest.root
