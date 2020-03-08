import multiprocessing as mp
import time



def testfunc(a):
    b = [2,3,4,a]
    time.sleep(30)
    return b


jobs = []
for i in range(4):
    p = mp.Process(target = testfunc, args=('a'))
    print('process',i,'started')
    jobs.append(p)
    time.sleep(0.1)
    p.start()

#wait for jobs to finish before continuing
for j in jobs:
    j.join()

print('all jobs done!')
print('jobs alive...')
print([j.is_alive() for j in jobs])