import taskdb
for i in range(10):
    T = taskdb.TaskDatabase('find_reps_%d' % i)
    T.clear_running()
    T.clear_runs()
    print T


