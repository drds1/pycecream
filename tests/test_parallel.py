import multiprocessing
import time
from random import randint

PROCESSES = 5
WORKER_CALLS = 7

def worker(num):
    """worker function"""
    print('Starting worker', num)
    time.sleep(randint(2,4))
    print('Exiting worker', num)
    return "ok"

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=PROCESSES)
    pool_outputs = pool.map(worker, range(WORKER_CALLS))
    pool.close()
    pool.join()
    print('Pool:', pool_outputs)