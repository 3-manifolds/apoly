import snappy, closed, taskdb, base64, sys


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
test_line15="K14n11024\tBdEDAAMDBAMEcrFOsRs5\t[[[-8, -3], [3, 1]]]"
test_line16="DT[obejfgHJImloKDCNaEb]\tBdEDAAIDBAQEOXIt5MZL\t[[[-5, -9], [-1, -2]]]"
test_line17="K14n15425\tBqYNAAEDAwQFBR6xkzlOTjk=\t[[[-5, 4], [1, -1]]]"
test_line18="DT[obfighjiMkaoelfcNDb]\tBtIOAAIEAwUEBUuTsbEbHjk=\t[[[-13, 8], [-5, 3]]]"
test_line19="DT[obcldEhjnLibgkoMFAc]\tBsQPAgMCBQQEBbE5Hk7kkzk=\t[[[-2, 7], [-1, 3]]]"
test_line20="K14n9577\tB2A/BAUCBAYGBgWNGznhTpyccg==\t[[[-3, 4], [-1, 1]]]"
test_line21="K14n17683\tB8I3AAUFBAMEBgYe0mPSHnI2OQ==\t[[[1, -3], [0, 1]]]"
test_line22="K14n23488\tB2A/BQIEBgYFBgVOHpOHNtKceA==\t[[[4, -7], [-1, 2]]]"
test_line23="DT[obghhdIkMjogCnabELf]\tB6Q+AQQEBgUFBgYbG5yxThuNjQ==\t[[[5, -4], [-1, 1]]]"
test_line24="DT[obcldfJbmgalNKCOeHI]\tB+A9AgUDBAYEBQZsGzk5sTlyHg==\t[[[-5, 3], [-2, 1]]]"
test_line25="DT[obcldEhknMJLbFGoIAc]\tCIT7AgYFBgcFBwcGLTYt5E6NbMZL\t[[[1, 4], [0, 1]]]"
test_line26="DT[obdkefkhaocMdLbnIGj]\tB8E+AAQFBAUGBgaNTrGNbHKT4Q==\t[[[8, -3], [3, -1]]]"
test_line27="DT[obdkeFGJLnHBKMCOIDa]\tCFHuAAIDBgYFBgcHjXLkGxvh5E5L\t[[[5, 3], [-2, -1]]]"
test_line28="DT[obghhcfIjbMgNlodEAk]\tBqQPAgMFBAQFBbFseGOHh3I=\t[[[5, -4], [-1, 1]]]"



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
    gluing_sols = closed.PHCGluingSolutionsOfClosed(M)
    return gluing_sols
    ans = gluing_sols.solutions()
    counts = map(len, ans)
    counts.append(sum(counts))
    new_data = [counts] + list(ans)
    return line + '\t' + '\t'.join(map(repr, new_data)), True

#find_reps(test_line1)

if __name__ == '__main__':
    n = int(sys.argv[1]) % 10
    T = taskdb.TaskDatabase('find_reps_%s' % n)
    print T
    T.run_function(find_reps, 50)
