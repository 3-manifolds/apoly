import taskdb, os
os.chdir('remaining')
db_names = ['find_reps_%d' % i for i in range(10)]
for db_name in db_names:
    T = taskdb.TaskDatabase(db_name)
    print T
    T.flat_file(db_name, where='done=0')
