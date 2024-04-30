import argparse
import multiprocessing
import time

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

class Generator:
  def __init__(
    self,
    TARGET: dict,
    CONSTRAINT: dict,
    PROCESS_LIMIT: int,
  ):
    self.TARGET = TARGET
    self.CONSTRAINT = CONSTRAINT

    self.ELEMENTS = list(TARGET.keys())
    self.LENGTH = sum(TARGET.values())

    self.PROCESS_LIMIT = PROCESS_LIMIT

  def generate(
    self,
    currentSequence: str,
    elements: list,
    length: int,
    strict: bool,
  ):
    sequences = []

    # Generation is complete if the sequence meets the target length
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
        sequences += self.generate(
          # Append the current element to the sequence
          currentSequence + element,
          # Filter available elements based on the target
          list(filter(
            lambda x: currentSequence.count(x) <= self.TARGET[x],
            elements
          )),
          length,
          strict,
        )

    return sequences

  def singleprocessGenerate(
    self,
  ):
    log('Generating sequences...', end='\r')

    sequences = self.generate(
      '',
      self.ELEMENTS,
      self.LENGTH,
      True,
    )

    log('Generated', len(sequences), 'sequences')

    return sequences

  def multiprocessGenerate(
    self,
  ):
    sequences = []

    # Find minimum length that fits the process limit
    startingLength = 0
    while len(self.ELEMENTS) ** startingLength < self.PROCESS_LIMIT:
      startingLength += 1

    log(
      'Using', len(self.ELEMENTS) ** startingLength, 'starting sequences',
      'and', self.PROCESS_LIMIT, 'processes'
    )

    # Generate starting sequences
    log('Generating sequences...', end='\r')
    startingSequences = self.generate(
      '',
      self.ELEMENTS,
      startingLength,
      False,
    )

    # Create a pool of processes
    pool = multiprocessing.Pool(processes=PROCESS_LIMIT)

    # Generate sequences in parallel
    results = pool.starmap(
      self.generate,
      ((
        startingSequence,
        self.ELEMENTS,
        self.LENGTH,
        True,
      ) for startingSequence in startingSequences)
    )

    for result in results:
      sequences += result

    log('Generated', len(sequences), 'sequences')

    return sequences

  def generateSequences(
    self,
  ):
    if (self.PROCESS_LIMIT > 1):
      return self.multiprocessGenerate()

    else:
      return self.singleprocessGenerate()

def main(
  TARGET: dict,
  CONSTRAINT: dict,
  OUTPUT_FILE: str,
  PROCESS_LIMIT: int,
):
  taskStart = time.time()

  generator = Generator(
    TARGET,
    CONSTRAINT,
    PROCESS_LIMIT,
  )

  sequences = generator.generateSequences()

  with open(OUTPUT_FILE, 'w') as file:
    for sequence in sequences:
      file.write(sequence + '\n')

  log('Done!', '(' + formatDuration(taskStart) + ')')

if __name__ == '__main__':
  parser=argparse.ArgumentParser()

  parser.add_argument('-t', '--target')
  parser.add_argument('-c', '--constraint')

  parser.add_argument('-p', '--processes')
  parser.add_argument('-o', '--output')
  parser.add_argument('-v', '--verbose', action='store_true')

  args=parser.parse_args()

  VERBOSE=args.verbose if args.verbose else False

  PROCESS_LIMIT=int(args.processes) if args.processes else 8
  OUTPUT_FILE=args.output if args.output else 'output.txt'

  CONSTRAINT = {}
  if args.constraint:
    for constraint in args.constraint.split(','):
      key, value = constraint.split('=')
      CONSTRAINT[key] = int(value)

  TARGET = {}
  if args.target:
    for target in args.target.split(','):
      key, value = target.split('=')
      TARGET[key] = int(value)

  main(
    TARGET,
    CONSTRAINT,
    OUTPUT_FILE,
    PROCESS_LIMIT,
  )
