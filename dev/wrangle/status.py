import taskdb
done = 0
for i in range(10):
    T = taskdb.TaskDatabase('find_reps_%d' % i)
    print T
    done += T.num_done_tasks()
print "DONE", done

