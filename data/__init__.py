import csv, collections

CuspedManifold = collections.namedtuple('CuspedManifold',
                                        ['name', 'L_space_fillings', 'finite_fillings', 'non_L_space_fillings',
                                         'thurston_norm', 'census_knot', 'longitude', 'id'])


def load_csv():
    file = open(__path__[0] + '/' + 'census_manifolds.csv')
    return csv.DictReader(file)

    
def csv_to_dataset():
    import dataset
    db = dataset.connect('sqlite:///:memory:')
    table = db['cusped']
    for data_dict in load_csv():
        table.insert(data_dict)
    table.create_index('id')
    table.create_index('names')
    return db, table
        
    

def eval_if_possible(item):
    try:
        return eval(item)
    except:
        return item
    
def load_data():
    ans = []
    for data_dict in load_csv():
        for key, value in data_dict.iteritems():
            data_dict[key] = eval_if_possible(value)
        ans.append(CuspedManifold(**data_dict))
    return ans

cusped_manifolds = load_data()
cusped_manifold_dict = collections.OrderedDict( [(M.name, M) for M in cusped_manifolds] )
