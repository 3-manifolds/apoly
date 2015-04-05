import taskdb, os
os.chdir('done')
done = 0
db_names = ['find_reps_%d' % i for i in range(10)]
for db_name in db_names:
    T = taskdb.TaskDatabase(db_name)
    print T
    T.flat_file(db_name, where='done=1')
    done += T.num_done_tasks()

print 'Done', done
os.system('cat ' + ' '.join(db_names) + '> find_reps_done_new')





