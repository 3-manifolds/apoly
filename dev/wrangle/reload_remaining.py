import taskdb

db_names = ['find_reps_%d' % i for i in range(10)]
conn = taskdb.connect_MySQL(None)
cur = conn.cursor()
for db in db_names:
    cur.execute('drop database ' + db)
cur.close(), conn.commit(), conn.close()

for db in db_names:
    taskdb.create_task_db('remaining/'+db, db)
    print taskdb.TaskDatabase(db)


