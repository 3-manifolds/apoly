import sys
sys.path.append('..')
import snappy
from manifold_reps import PSL2RRepOf3ManifoldGroup
from manifold_reps.real_reps import shift_of_central

sample_data = [
               ('m004(3,2)', [0.48886560625734599, 0.25766090533555303]), 
               ('m030(1,3)', [-3.0413531153050002, -0.39791078825716603, -1.0980127709735299, -1.5128763968640999]),
               ('m032(5,2)', [1.4397193529227299, -0.37244727986758702, 0.67639598574213, 0.249420323040287]),
               ('m032(7,1)', [-0.53149787942290505, 1.88940106877921, -0.367853983707157, 0.87494072873215301]),
               ('m034(-2,3)', [0.64486997829249104, 5.2944764441563397, 0.255973524222671, -0.044117577967982401]),
               ('m071(-4, 1)', [(0.55948912617811042+0.23400162016416198j), (0.67650039034242671-1.1474257693865324e-17j), (1.7704952148495239-0.94049607714155237j), (2.643772147889158+1.5020248611456067e-15j), (0.2114871099765474-2.367202429557164e-16j)]),
               ('m071(-5,1)', [(0.55276546381527392+0.17809590030216793j), (0.68710050225984587-2.043928395675953e-17j), (1.9299238691198983-0.76852635737322539j), (2.5966309123403892+1.6688390551847239e-15j), (0.21585359594747663-2.7103458327357865e-16j)]),
                ('m071(-6,1)', [(0.54967530746810289+0.14417261050664634j), (0.6923450301230577-9.3067964391777209e-20j), (2.0141720155237239-0.64484235997545269j), (2.5738125283754054+3.8720715076017058e-16j), (0.21803368751082777-6.5573140253643421e-17j)]),
                 ('m071(-7,1)', [(0.54798542957115404+0.12127149741227804j), (0.6953255832254438-1.271712523437961e-17j), (2.0637677963047372-0.5536905828758234j), (2.5609893276405575+1.059570683450395e-15j), (0.21927860212419115-1.7743593858656029e-16j)]),
                  ('m071(-8,1)', [(0.54695606619291637+0.10472467858143307j), (0.697182863730199+1.8235651156127249e-17j), (2.0953296664174421-0.4843519788329127j), (2.5530510662630088-4.5116226623302756e-16j), (0.22005654637331931+7.2018716811597148e-17j)]),
                ('m071(-9,1)', [(0.54628089812013647+0.092189810267817446j), (0.69841893915573283+5.5276413020141601e-18j), (2.1166221992371859-0.43007005468329118j), (2.5477899694764905-1.1178170508975334e-16j), (0.22057523232661366+1.760928580844308e-17j)]),
                 ('m071(-10,1)', [(0.54581340342992435+0.082356146274548167j), (0.69928333205054738+4.382286003266566e-18j), (2.1316509189163386-0.3865251775598425j), (2.5441212712396375-2.6221343423297572e-17j), (0.22093840030687181+3.1696831381799638e-18j)]),
                ('m071(-31,1)',[(0.54415721128853944+0.025536087441187939j), (0.70240372066354295-1.4011788256507199e-18j), (2.186875951813898-0.12250770860367634j), (2.530948300574003-1.6472990780420066e-17j), (0.22225249686975052+3.3364726965139508e-18j)]),
                ('m071(-32,1)', [(0.54414627241049263+0.024725258122055737j), (0.70242463990219672-1.2840841093818995e-18j), (2.1872514154953113-0.11863532653757804j), (2.5308603599283024+3.709569138633798e-17j), (0.22226132303461835-6.1092113139711614e-18j)]),
    ('m071(-50,1)', [(0.54976877948338143-0.14532673525108694j), (0.69218257003081918-1.4880788668469433e-17j), (2.0115052476293531+0.64927858677364403j), (2.5745144729398661+2.6529452443703887e-19j), (0.21796595553894438+4.6409889096230281e-18j)]),
    ('m071(1,7)', [1.372050472100670 + 0.8140025989721413j, 0.3147281566257106 - 9.146398497790777e-18j, -0.464469854756940 - 1.01620531962142j, -1.159134584846391 - 2.608331729328281e-17j, 2.364812291169813 - 3.770349594140284e-17j])
]

def sample_rep(index, precision):
    name, shapes = sample_data[index]
    M = snappy.Manifold(name)
    return PSL2RRepOf3ManifoldGroup(M, shapes, precision, [True, False, True])

def sample_rep_with_peripheral_renormalization(index, precision):
    name, shapes = sample_data[index]
    M = snappy.Manifold(name)
    M.set_peripheral_curves('fillings')
    return PSL2RRepOf3ManifoldGroup(M, 
                                    target_meridian_holonomy_arg=0.0, 
                                    rough_shapes=shapes,
                                    precision=precision,
                                    fundamental_group_args=[True, False, True])
test_bits = 256

def basic_test():
    for i in range(len(sample_data)):
        rho = sample_rep(i, test_bits)
        error = rho.polished_holonomy().check_representation()
        print rho.manifold, rho.manifold.homology(), error.log(2).ceil(), rho.representation_lifts()

def current_test():
    for i in [k for k in range(len(sample_data)) if not k in [1]]:
        rho = sample_rep_with_peripheral_renormalization(i, test_bits)
        meridian, longitude = rho.polished_holonomy().peripheral_curves()[0]
        error = rho.polished_holonomy().check_representation()
        print sample_data[i][0],  rho.manifold.homology(), error.log(2).ceil(), rho.representation_lifts()
        print rho.manifold
        print rho.manifold.cusp_info()
        rho_til = rho.lift_on_cusped_manifold()
        print "   euler cocycle: ", [-shift_of_central(rho_til(R)) for R in rho.relators()]
        print "   cobdr 1: ", repr(rho.coboundary_1_matrix()).replace('\n', '\n' + 13*' ')
        print "   new meridian: ", rho_til(meridian)
        print "   new longitude: ", rho_til(longitude)

        # rho =  sample_rep(i, test_bits)
        # meridian, longitude = rho.polished_holonomy().peripheral_curves()[0]
        # if abs(rho(meridian).trace()) <= 2 and abs(rho(longitude).trace()) <=2:
        #     p, q = map(int, rho.manifold.cusp_info(0).filling)
        #     rho_til = rho.lift_on_cusped_manifold()
        #     m_shift = translation_amount(rho_til(meridian))
        #     l_shift = translation_amount(rho_til(longitude))
        #     combo = p*m_shift + q*l_shift
        #     print "   old meridian trans: ", RDF(m_shift)
        #     print "   old longitude trans: ", RDF(l_shift)
        #     print "   shift of filling curve: ", RDF(combo)
        print "\n"


    

current_test()

#rho = sample_rep(-1, 1000)
#D = lift_on_cusped_manifold(rho)

#print rho.polished_holonomy()



  
