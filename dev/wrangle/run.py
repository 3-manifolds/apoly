from lift_test import save_data
import taskdb

def run_line(line):
    try: 
        save_data(line)
        return line, True
    except:
        return line, False

TD = taskdb.TaskDatabase("SLarcs2")
TD.run_function(run_line, num_tasks=1)

