import argparse
import multiprocessing
import time
import math
from tqdm import tqdm
import random

PROCESS = 0
VERBOSE = False

def log(*args, end='\n'):
  if VERBOSE:
    print(*args, end=end)

def formatDuration(
    startTime: float,
):
  duration = time.time() - startTime

  if duration < 1:
    return f'{duration * 1000:.2f}ms'
  elif duration < 60:
    return f'{duration:.2f}s'
  else:
    return f'{duration // 60:.0f}m {duration % 60:.2f}s'

def taskWrapper(
    pool: multiprocessing.Pool,
    handler: callable,
    tasks: list,
    unit = 'tasks'
):
  iterator = pool.imap(
    handler,
    (task for task in tasks)
  )

  if VERBOSE:
    return list(tqdm(
      iterator,
      total=len(tasks),
      unit=' ' + unit,
      bar_format='{percentage:3.0f}% | {n_fmt}/{total_fmt} {unit} | {rate_fmt} | {elapsed}<{remaining}'
    ))
  else:
    return list(iterator)

class Filter:
  def __init__(
    self,
    INPUT_FILE: str,
    OUTPUT_FILE: str,

    DISTANCE: float,
    ITERATION: int,
    CHUNK: int,
  ):
    self.INPUT_FILE = INPUT_FILE
    self.OUTPUT_FILE = OUTPUT_FILE

    self.DISTANCE = DISTANCE
    self.ITERATION = ITERATION
    self.CHUNK = CHUNK

  def __distance__(
      self,
      str1: str,
      str2: str,
  ):
    distance = 0
    for i in range(len(str1)):
      if str1[i] != str2[i]:
        distance += 1
    return distance

  def __filter__(
    self,
    params: list
  ):
    sequences = params[0]
    currentFiltered = params[1]

    newFiltered = []

    filtered = list(currentFiltered) \
      if len(currentFiltered) \
      else newFiltered

    for sequence in sequences:
      if all(self.__distance__(filtered[i], sequence) >= self.DISTANCE
        for i in range(len(filtered))
      ):
        newFiltered.append(sequence)

    return newFiltered

  def run(self):
    sequences = []
    with open(INPUT_FILE, 'r') as file:
      sequences = file.read().splitlines()

    pool = multiprocessing.Pool(processes=PROCESS)

    for i in range(self.ITERATION):
      log('Filtering', len(sequences), 'sequences...')

      chunkCount = math.ceil(len(sequences) / self.CHUNK)

      # Create tasks
      tasks = []
      for i in range(0, chunkCount):
        tasks.append([
          sequences[i * self.CHUNK:(i + 1) * self.CHUNK],
          []
        ])

      # Execute tasks
      results = taskWrapper(
        pool,
        self.__filter__,
        tasks,
        'chunks'
      )

      # Merge results
      sequences = []
      for result in results:
        sequences += result

      # Write sequences to file
      with open(OUTPUT_FILE, 'w') as file:
        for sequence in sequences:
          file.write(sequence + '\n')

      if (len(results[0]) == self.CHUNK or chunkCount == 1):
        break

    log('Filtered', len(sequences), 'sequences')

    # Write sequences to file
    with open(OUTPUT_FILE, 'w') as file:
      sequences.sort()
      for sequence in sequences:
        file.write(sequence + '\n')

    return sequences


def main(
  INPUT_FILE: str,
  OUTPUT_FILE: str,
  DISTANCE: str,
  ITERATION: int,
  CHUNK: int,
):
  log('Starting...', end='\r')
  taskStart = time.time()

  filter = Filter(
    INPUT_FILE,
    OUTPUT_FILE,
    DISTANCE,
    ITERATION,
    CHUNK,
  )

  filter.run()

  log('Done!', '(' + formatDuration(taskStart) + ')')

if __name__ == '__main__':
  parser=argparse.ArgumentParser(
    epilog='> python3 filter.py --input input.txt --output output.txt --distance 4 --iteration 5 --chunk 1000 --process 8 --verbose'
  )

  parser.add_argument('-i', '--input', required=True)
  parser.add_argument('-o', '--output', required=True)

  parser.add_argument('-d', '--distance', required=True, help='Minimum distance between sequences')
  parser.add_argument('-n', '--iteration', required=True, help='Number of iterations to run')
  parser.add_argument('-c', '--chunk', required=True, help='Size of chunks to process')

  parser.add_argument('-p', '--process')
  parser.add_argument('-v', '--verbose', action='store_true')

  args=parser.parse_args()

  INPUT_FILE=args.input
  OUTPUT_FILE=args.output

  DISTANCE = float(args.distance)
  ITERATION = int(args.iteration) if args.iteration else 0
  CHUNK = int(args.chunk) if args.chunk else 0

  PROCESS=int(args.process) if args.process else 1
  VERBOSE=args.verbose if args.verbose else False

  main(
    INPUT_FILE,
    OUTPUT_FILE,
    DISTANCE,
    ITERATION,
    CHUNK,
  )
