from subprocess import Popen
from subprocess import DEVNULL
from multiprocessing import Pool
import multiprocessing
import glob
import os
import sys
import getopt


def worker(command):
   p = Popen(command, stdout=DEVNULL, stderr=DEVNULL, shell=True)
   p.wait()

def main(argv):
   maxProcessors = None
   trialsPerIstance = 10
   exePath = ""
   instancesPath = ""

   try:
      opts, args = getopt.getopt(argv, "e:i:p:t:", ["exe=", "ipath=", "procs=", "trials="])
   except getopt.GetoptError:
      print("doe.py -e <exe_PATH> -i <input_PATH> [-p <N_PROCESSORS>] [-t <N_TRIALS>]")
      sys.exit(-1)

   required = ["-e", "-i"]
   optStrings = [opt[0] for opt in opts]

   for req in required:
      if req not in optStrings:
         print("doe.py -e <exe_PATH> -i <input_PATH> [-p <N_PROCESSORS>] [-t <N_TRIALS>]")
         sys.exit(-1)

   for opt, arg in opts:
      if opt == '-e':
         exePath = arg
      elif opt == "-i":
         instancesPath = arg
      elif opt == "-p":
         maxProcessors = int(arg)
      elif opt == "-t":
         trialsPerIstance = int(arg)

   maxProcessors = maxProcessors if maxProcessors is not None else multiprocessing.cpu_count()

   print("Il test verrà eseguito con i seguenti parametri:")
   print("Percorso eseguibile: {0}".format(exePath))
   print("Percorso istanze: {0}".format(instancesPath))
   print("Processori utilizzabili: {0}".format(str(maxProcessors)))
   print("Tentativi per istanza: {0}".format(str(trialsPerIstance)))

   nfesArray = [100000]
   hArray = [100]
   tArray = [10]
   alphaMinArray = [0.2, 0.25, 0.3, 0.35, 0.4]
   alphaMaxArray = [0.6, 0.65, 0.7, 0.75, 0.8]

   args = []

   #create exe params
   for nfes in nfesArray:
      for h in hArray:
         for t in tArray:
            for alphaMin in alphaMinArray:
               for alphaMax in alphaMaxArray:
                  args.append("{0} {1} {2} {3} {4}".format(nfes, h, t, alphaMin, alphaMax))

   instances = glob.glob(os.path.join(instancesPath, "*.txt"))
   for index, instance in enumerate(instances):
      instances[index] = instance.replace("\\", "/")

   pool = Pool(maxProcessors)
   results = []

   for index, funcParams in enumerate(args):
      commands = ["{0} {1} {2}".format(exePath, instance, funcParams) for instance in instances]
      for command in commands:
         for i in range(trialsPerIstance):
            results.append(pool.apply_async(worker, args=(command,)))

   pool.close()
   pool.join()

if __name__ == '__main__':
    main(sys.argv[1:])