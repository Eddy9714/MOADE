from subprocess import Popen
from subprocess import DEVNULL, STDOUT
import glob
import os
import sys
import getopt
import math
import csv
import shutil

def main(argv):
   maxProcessors = 30
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

   print("Il test verrà eseguito con i seguenti parametri:")
   print("Percorso eseguibile: {0}".format(exePath))
   print("Percorso istanze: {0}".format(instancesPath))
   print("Processori utilizzabili: {0}".format(str(maxProcessors)))
   print("Tentativi per istanza: {0}".format(str(trialsPerIstance)))

   nfesArray = [100000]
   hArray = [49, 99, 149]
   tArray = [5, 10, 15, 20]
   alphaMinArray = [0., 0.1, 0.2, 0.3, 0.4]
   alphaMaxArray = [0.6, 0.7, 0.8, 0.9, 1.]

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

   processorsAvailable = maxProcessors
   procs = []

   for index, funcParams in enumerate(args):
      print("Eseguo con i parametri: {0}".format(funcParams))

      commands = ["{0} {1} {2}".format(exePath, instance, funcParams) for instance in instances]

      for command in commands:
         trialsRemaining = trialsPerIstance

         for i in range(trialsRemaining):
            procs.append(Popen(command, stdout=DEVNULL, stderr=DEVNULL, shell=True))
            processorsAvailable = processorsAvailable - 1

            if processorsAvailable == 0:
               procs[0].wait()
               processorsAvailable = processorsAvailable + 1
               procs = procs[1:]

      for proc in procs:
         proc.wait()

      if index > 0:

         folders = [directory for directory in os.listdir(instancesPath) if os.path.isdir(os.path.join(instancesPath, directory))]
         bucket = dict()
         score = [0., 0.]
         worstKey = None

         for folder in folders:
            keys = folder.split("-")
            mainKey = keys[0]
            if bucket.get(mainKey) is None:
               bucket[mainKey] = []

            bucket[mainKey].append(folder)

         for key, value in bucket.items():
            first_path = os.path.join(instancesPath, value[0])
            second_path = os.path.join(instancesPath, value[1])

            x1, y1 = paretoFrontFromPath(first_path)
            x2, y2 = paretoFrontFromPath(second_path)
            x3, y3 = optimalParetoFrontFromValues(x1, y1, x2, y2)

            max_x3 = max(x3)
            max_y3 = max(y3)

            refPoint = (1.1 * max_x3, 1.1 * max_y3)
            hv_moade1 = hypervolume(x1, y1, refPoint)
            hv_moade2 = hypervolume(x2, y2, refPoint)

            if hv_moade1 > hv_moade2:
               score[0] += 1
            elif hv_moade2 > hv_moade1:
               score[1] += 1

            igd_moade1 = igd(x1, y1, x3, y3)
            igd_moade2 = igd(x2, y2, x3, y3)

            if igd_moade1 < igd_moade2:
               score[0] += 1
            elif igd_moade2 < igd_moade1:
               score[1] += 1

            gd_moade1 = gd(x1, y1, x3, y3)
            gd_moade2 = gd(x2, y2, x3, y3)

            if gd_moade1 < gd_moade2:
               score[0] += 1
            elif gd_moade2 < gd_moade1:
               score[1] += 1

            cov_moade1 = cov_two_sets(x1, y1, x2, y2)
            cov_moade2 = cov_two_sets(x2, y2, x1, y1)

            if cov_moade1 > cov_moade2:
               score[0] += 1
            elif cov_moade2 > cov_moade1:
               score[1] += 1

            if score[0] < score[1]:
               worstKey = value[0].split("-")[1]
            else:
               worstKey = value[1].split("-")[1]

         for key, value in bucket.items():
            for el in value:
               if worstKey in el:
                  shutil.rmtree(os.path.join(instancesPath, el), ignore_errors=True)


def cov_two_sets(x1, y1, x2, y2):
   count = 0

   for i in range(len(x2)):
      for j in range(len(x1)):
         if x1[j] <= x2[i] and y1[j] <= y2[i] and (x1[j] < x2[i] or y1[j] < y2[i]):
            count += 1
            break

   return count / len(x2)


def gd(x1, y1, x_best, y_best):
   v = 0
   for i in range(len(x1)):
      min_distance = math.inf
      for j in range(len(x_best)):
         distance = math.sqrt(pow(x_best[j] - x1[i], 2) + pow(y_best[j] - y1[i], 2))
         if distance < min_distance:
            min_distance = distance

      v += pow(min_distance, 2)

   v = math.sqrt(v) / len(x1)

   return v


def igd(x1, y1, x_best, y_best):
   v = 0
   for i in range(len(x_best)):
      min_distance = math.inf
      for j in range(len(x1)):
         distance = math.sqrt(pow(x_best[i] - x1[j], 2) + pow(y_best[i] - y1[j], 2))
         if distance < min_distance:
            min_distance = distance

      v += pow(min_distance, 2)

   v = math.sqrt(v) / len(x_best)

   return v


def hypervolume(x, y, ref):
   area = 0

   for i in range(1, len(x)):
      w = x[i] - x[i - 1]
      h = ref[1] - y[i - 1]
      area += w * h

   area += (ref[0] - x[len(x) - 1]) * (ref[1] - y[len(x) - 1])

   return area


def optimalParetoFrontFromValues(x1, y1, x2, y2):
   x = x1 + x2
   y = y1 + y2

   points = list(zip(x, y))

   return paretoFrontFromPoints(points)


def paretoFrontFromPath(path):
   pts = loadDataFromPath(path)
   return paretoFrontFromPoints(pts)


def loadDataFromPath(path):
   results = glob.glob(os.path.join(path, "*.csv"))
   pts = []

   for result in results:
      with open(result) as file:
         csv_reader = csv.reader(file, delimiter=";")
         for i in range(3):
            next(csv_reader)
         for row in csv_reader:
            pts.append([float(row[0]), float(row[1])])

   return pts


def avgTimesFromPath(path):
   results = glob.glob(os.path.join(path, "*.csv"))
   tot = 0
   for result in results:
      with open(result) as file:
         csv_reader = csv.reader(file, delimiter=";")
         csvl = list(csv_reader)
         tot += float(csvl[1][5])

   return tot / len(results)


def paretoFrontFromPoints(points):
   result = findParetoDominantPoints(points)

   x = []
   y = []
   for point in result:
      x.append(point[0])
      y.append(point[1])

   return x, y


def findParetoDominantPoints(pts):
   result = []
   pts = sorted(pts, key=lambda element: element[0])

   i = 0
   j = 0

   while i < len(pts):
      result.append(pts[i])

      j = i + 1
      while j < len(pts):
         if pts[j][1] < pts[i][1]:
            break
         j += 1
      i = j

   return result


if __name__ == '__main__':
    main(sys.argv[1:])