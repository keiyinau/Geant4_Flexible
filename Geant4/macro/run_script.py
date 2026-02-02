import os 
import glob

filenames=glob.glob("SquareBox*.txt")
for filename in filenames:
    file_run=filename.split(".")[0]
    with open("run_script_%s.mac"%(file_run), "w") as file:
        file.write("/run/numberOfThreads 20\n")
        file.write("/run/initialize\n")
        file.write("/control/verbose 0\n")
        file.write("#/gps/verbose 2\n")
        file.write("/gps/particle o-Ps\n")
        file.write("/gps/pos/type Volume\n")
        file.write("/gps/pos/shape Cylinder\n")
        file.write("/gps/pos/centre 0 0 0 mm\n")
        file.write("/gps/pos/radius 6 mm\n  ")
        file.write("/gps/pos/halfz 2 mm\n")
        file.write("/gps/hist/type biaspt\n")
        file.write("/gps/hist/file biaspt.dat\n")
        file.write("/gps/hist/type biaspp\n  ")
        file.write("/gps/hist/file biaspp.dat\n")
        file.write("/gps/hist/type biasz\n")
        file.write("/gps/hist/file biasz.dat\n")
        file.write("/gps/ene/type Mono\n")
        file.write("/gps/ene/mono 0. MeV\n")
        file.write("/gps/number 1\n")
        file.write("/gps/ang/type iso\n")
        file.write("/MyDetector/setFileName Coordinate_file_%s\n"%(filename.split(".")[0]))
        file.write("/MyDetector/setDetectorCoordinate %s\n"%(filename))
        file.write("/run/beamOn 100000")
    os.system("./sim run_script_%s.mac >> output_%s.log"%(file_run,file_run))