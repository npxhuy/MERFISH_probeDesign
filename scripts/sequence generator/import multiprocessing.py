import multiprocessing
import time

def generate(values, length, current=""):
    results = []

    if len(current) == length:
        if 'GGGG' not in current:
            results.append(current)
    else:
        for value in values:
            results += generate(values, length, current + value)

    return results

def splitGeneration(values, length, limit):
    results = []
    baseValues = []

    if limit:
        # If limit is specified, generate base values
        baseValues = generate(values, limit)
        # Create a pool of worker processes
        with multiprocessing.Pool() as pool:
            # Use starmap to apply generate function to each base value
            results = pool.starmap(generate, [(values, length, baseValue) for baseValue in baseValues])
    else:
        # If no limit is specified, generate all combinations
        results = generate(values, length)

    print(f"Generated {len(results)} combinations from {len(baseValues)} base values.")
    return results

if __name__ == '__main__':
    start_time = time.time()
    results = splitGeneration(["G"]*5 + ["A"]*3 + ["T"]*3, 11, 0)
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds") 