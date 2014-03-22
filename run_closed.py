import snappy, closed, taskdb, base64


test_line = "DT[obejfgjMkeaonldbCih]\tAg4BAQF4G9I=\t[[[21, -4], [-5, 1]]]"

def manifold_from_bytes_n_cobs(encoded_bytes, cobs):
    R = snappy.Manifold('empty')
    R._from_bytes(base64.decodestring(encoded_bytes))
    R.set_peripheral_curves('combinatorial')
    R.set_peripheral_curves(cobs)
    R.dehn_fill( (1,0) )
    return R

def find_reps(line):
    line = line.strip()
    name, encoded_bytes, cobs = line.split('\t')
    M = manifold_from_bytes_n_cobs(encoded_bytes, eval(cobs))
    ans = closed.PHCGluingSolutionsOfClosed(M).solutions()
    counts = map(len, ans)
    counts.append(sum(counts))
    new_data = [counts] + list(ans)
    return line + '\t' + '\t'.join(map(repr, new_data)), True

T = taskdb.TaskDatabase('find_reps')
T.run_function(find_reps, 50)
