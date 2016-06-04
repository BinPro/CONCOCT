"""
Processes and parallelism
=========================

The ``processes`` module provides some convenience functions
for using processes in python.

Adapted from http://stackoverflow.com/a/16071616/287297

Example usage:

    print prll_map(lambda i: i * 2, [1, 2, 3, 4, 6, 7, 8], 32, verbose=True)

Comments:

"It spawns a predefined amount of workers and only iterates through the input list
 if there exists an idle worker. I also enabled the "daemon" mode for the workers so
 that KeyboardInterupt works as expected."

Pitfalls: all the stdouts are sent back to the parent stdout, intertwined.

Alternatively, use this fork of multiprocessing: https://github.com/uqfoundation/multiprocess
"""

# Modules #
import multiprocessing
from tqdm import tqdm

################################################################################
def target_func(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None: break
        q_out.put((i, f(x)))

################################################################################
def prll_map(func_to_apply, items, cpus=None, verbose=False):
    # Number of cores #
    if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
    # Create queues #
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()
    # Process list #
    new_proc = lambda t,a: multiprocessing.Process(target=t, args=a)
    processes = [new_proc(target_func, (func_to_apply, q_in, q_out)) for x in range(cpus)]
    # Start them all #
    for proc in processes:
        proc.daemon = True
        proc.start()
    # Put and get #
    sent = [q_in.put((i, x)) for i, x in enumerate(items)]
    for x in range(cpus): q_in.put((None, None))
    res  = [q_out.get() for x in range(len(sent))]
    # Wait for them to finish #
    if verbose:
        for proc in tqdm(processes): proc.join()
    else:
        for proc in processes: proc.join()
    # Return results #
    return [x for i, x in sorted(res)]