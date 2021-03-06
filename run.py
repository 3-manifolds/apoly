#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --nice=10000
#SBATCH --time=1-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_out/%j

import slurm
import taskdb2
import snappy
import census_data
import lift_test
import traceback


def create_database():
    names = [M.name for M in census_data.cusped_manifolds]
    taskdb2.create_task_db(names, 'PEChar', overwrite=True)
    db = taskdb2.TaskDatabase('PEChar')
    db.add_column('lifter', 'mediumblob')
    return db

def add_lifter(task):
    try:
        task['lifter'] = lift_test.save_data(task['name'], to_file=False)
        task['done'] = True
    except:
        print('<apoly-error>')
        print('task: ' + task['name'])
        print(traceback.format_exc())
        print('</apoly-error>')


if __name__ == '__main__':
    db = taskdb2.TaskDatabase('PEChar')
    db.run_function(add_lifter, num_tasks=1)

