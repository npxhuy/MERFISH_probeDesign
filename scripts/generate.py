import argparse
import multiprocessing
import time
from tqdm import tqdm

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

class Generator:
  def __init__(
    self,
    OUTPUT_FILE: str,

    TARGET: dict,
    CONSTRAINT: dict,
  ):
    self.OUTPUT_FILE = OUTPUT_FILE

    self.TARGET = TARGET
    self.CONSTRAINT = CONSTRAINT

    self.ELEMENTS = list(TARGET.keys())
    self.LENGTH = sum(TARGET.values())

  def __generate__(
    self,
    params: list,
  ):
    # Parse parameters
    currentSequence = params[0]
    elements = params[1]
    length = params[2]
    strict = params[3]

    sequences = []

    # Generation is complete if the sequence meets the target count
    if len(currentSequence) == length:

      # Strict verification checks if the sequence meets the target
      if strict:
        if all(currentSequence.count(element) == self.TARGET[element]
            for element in self.TARGET.keys()) \
          and all(currentSequence.count(subsequence) == self.CONSTRAINT[subsequence]
            for subsequence in self.CONSTRAINT.keys()
        ):
          sequences.append(currentSequence)

      # Non-strict verification checks if the sequence is within the target
      else:
        if all(currentSequence.count(element) <= self.TARGET[element]
            for element in self.TARGET.keys()) \
          and all(currentSequence.count(subsequence) <= self.CONSTRAINT[subsequence]
            for subsequence in self.CONSTRAINT.keys()
        ):
          sequences.append(currentSequence)

    # Keep generating if the sequence is not complete
    else:
      for element in elements:
        sequences += self.__generate__([
          # Append the current element to the sequence
          currentSequence + element,
          # Filter available elements based on the target
          list(filter(
            lambda x: currentSequence.count(x) <= self.TARGET[x],
            elements
          )),
          length,
          strict,
        ])

    return sequences

  def run(self):
    # Find minimum length that fits the process limit
    startingLength = 0
    while len(self.ELEMENTS) ** startingLength < PROCESS:
      startingLength += 1

    # Generate starting sequences
    log('Generating initial sequences...')
    startingSequences = self.__generate__([
      '',
      self.ELEMENTS,
      startingLength,
      False,
    ])

    # Generate sequences in parallel
    log('Generating final sequences...')
    pool = multiprocessing.Pool(processes=PROCESS)

    # Create tasks
    tasks = []
    for startingSequence in startingSequences:
      tasks.append([
        startingSequence,
        self.ELEMENTS,
        self.LENGTH,
        True,
      ])

    # Execute tasks
    results = taskWrapper(
      pool,
      self.__generate__,
      tasks,
      'branches'
    )

    # Merge results
    sequences = []
    for result in results:
      sequences += result

    log('Generated', len(sequences), 'sequences')

    # Write sequences to file
    with open(self.OUTPUT_FILE, 'w') as file:
      for sequence in sequences:
        file.write(sequence + '\n')

    return sequences

def main(
  TARGET: dict,
  CONSTRAINT: dict,
  OUTPUT_FILE: str,
):
  log('Starting...', end='\r')
  taskStart = time.time()

  generator = Generator(
    OUTPUT_FILE,
    TARGET,
    CONSTRAINT,
  )

  generator.run()

  log('Done!', '(' + formatDuration(taskStart) + ')')

if __name__ == '__main__':
  parser=argparse.ArgumentParser(
    epilog='> python3 generate.py --output output.txt --target A=5,G=10,T=5 --constraint GGGG=0 --process 8 --verbose'
  )

  parser.add_argument('-o', '--output', required=True)

  parser.add_argument('-t', '--target', required=True, help='Target count for each element. Example: A=5,G=10,T=5')
  parser.add_argument('-c', '--constraint', required=False, help='Constraint count for subsequences. Example: GGGG=0')

  parser.add_argument('-p', '--process')
  parser.add_argument('-v', '--verbose', action='store_true')

  args=parser.parse_args()

  OUTPUT_FILE=args.output

  TARGET = {}
  for target in args.target.split(','):
    key, value = target.split('=')
    TARGET[key] = int(value)

  CONSTRAINT = {}
  if args.constraint:
    for constraint in args.constraint.split(','):
      key, value = constraint.split('=')
      CONSTRAINT[key] = int(value)

  PROCESS=int(args.process) if args.process else 1
  VERBOSE=args.verbose if args.verbose else False

  main(
    TARGET,
    CONSTRAINT,
    OUTPUT_FILE
  )
