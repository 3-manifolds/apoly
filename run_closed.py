import snappy, closed, taskdb, base64


test_line0 = "DT[obejfgjMkeaonldbCih]\tAg4BAQF4G9I=\t[[[21, -4], [-5, 1]]]"
test_line1 = "L13n5884\tAzoBAQICNpwejQ==\t[[[5, -1], [1, 0]]]"
test_line2 = "L14n63017\tBOwBAgMDAxuxGxtO\t[[[3, 4], [-1, -1]]]"
test_line3 = "L13n3710\tBPgCAQIDA4dLHuGN\t[[[7, -4], [2, -1]]]"
test_line4 = "DT[occdhdEHiLkofCmnBajg]\tBOoBAgMDAxtOGxux\t[[[-4, -11], [-1, -3]]]"
test_line5 = "L14n53729\tBOwCAQMDA7GxG06x\t[[[-4, -3], [-1, -1]]]"
test_line6 = "L14n55098\tBPEAAwIDA41OjbEb\t[[[11, -4], [3, -1]]]"
test_line7 = "DT[obdkefhlijObnacmdkG]\tBfADAgIDAwQEjWxybB45\t[[[7, 19], [-3, -8]]]"
test_line8="K14n26950\tCVHmAwACAwYGBwYIBwg5cuQbG07kTktL\t[[[3, -4], [1, -1]]]"
test_line9="L14n27443\tDMBt6wIFBAYHCAUHCQoLCwuTTpMteE6THuSxG06x\t[[[-1, 1], [-1, 0]]]"
test_line10="DT[obcldefcijLMKabONGH]\tBdEDAAIDBAMEjeEb5B45\t[[[-19, -8], [-7, -3]]]"
test_line11="K14n18445\tBvAOAQIEAwUEBTke4eFOOUs=\t[[[7, 5], [-3, -2]]]"
test_line12="DT[obejfDHLKaoICnEBgjm]\tBsIPAAQFAwMEBY0b5BtyjXI=\t[[[8, -3], [3, -1]]]"
test_line13="L14n14006\tBsgPAQQDAwUFBU6ckx5jcpM=\t[[[5, 8], [-2, -3]]]"
test_line14 = "L14n32374\tBsEPAAQDBQMEBY3kGxuNcnI=\t[[[-3, -1], [1, 0]]]"
test_line15= "K14n11024\tBdEDAAMDBAMEcrFOsRs5\t[[[-8, -3], [3, 1]]]"
test_line16= "DT[obejfgHJImloKDCNaEb]\tBdEDAAIDBAQEOXIt5MZL\t[[[-5, -9], [-1, -2]]]"

def manifold_from_bytes_n_cobs(encoded_bytes, cobs):
    R = snappy.Manifold('empty')
    R._from_bytes(base64.decodestring(encoded_bytes))
    R.set_peripheral_curves('combinatorial')
    R.set_peripheral_curves(cobs)
    R.dehn_fill( (1,0) )
    return R

def find_reps(line):
    line = line.strip()
    print line
    name, encoded_bytes, cobs = line.split('\t')
    M = manifold_from_bytes_n_cobs(encoded_bytes, eval(cobs))
    ans = closed.PHCGluingSolutionsOfClosed(M).solutions()
    counts = map(len, ans)
    counts.append(sum(counts))
    new_data = [counts] + list(ans)
    return line + '\t' + '\t'.join(map(repr, new_data)), True

#find_reps(test_line16)
T = taskdb.TaskDatabase('find_reps')
T.run_function(find_reps, 50)


