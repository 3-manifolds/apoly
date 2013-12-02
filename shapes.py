from snappy.snap.shapes import pari, gen, float_to_pari, complex_to_pari, pari_column_vector
from snappy.snap.shapes import enough_gluing_equations, eval_gluing_equation
from snappy.snap.shapes import infinity_norm, pari_matrix, pari_vector_to_list, _within_sage, ComplexField
from sage.all import exp, CC

def gluing_equation_errors(eqns, shapes, RHS_of_last_eqn):
    last = [eval_gluing_equation(eqns[-1], shapes) - RHS_of_last_eqn]
    return [eval_gluing_equation(eqn, shapes) - 1 for eqn in eqns[:-1]] + last

def gluing_equation_error(eqns, shapes, RHS_of_last_eqn):
    return infinity_norm(gluing_equation_errors(eqns, shapes, RHS_of_last_eqn))

def sage_complex_to_pari(z, dec_prec):
    return pari.complex( float_to_pari(z.real(), dec_prec), float_to_pari(z.imag(), dec_prec) )


def polished_tetrahedra_shapes(manifold, target_meridian_log_holonomy,
                               dec_prec=None, bits_prec=200, ignore_solution_type=False):
    """
    Refines the current solution to the gluing equations to one with
    the specified accuracy.  
    """

    if dec_prec is None:
        dec_prec = gen.prec_bits_to_dec(bits_prec)
    else:
        bits_prec = gen.prec_dec_to_bits(dec_prec)
    working_prec = dec_prec + 10
    target_espilon = float_to_pari(10.0, working_prec)**-dec_prec
    
    # This is a potentially long calculation, so we cache the result

    if "polished_shapes" in manifold._cache.keys():
        curr_sol = manifold._cache["polished_shapes"]
        if curr_sol[0].precision() >= gen.prec_dec_to_words(dec_prec):
            if _within_sage:
                CC = ComplexField(bits_prec)
                return [CC(z) for z in curr_sol]
            return [s.precision(dec_prec) for s in curr_sol]

    init_shapes = pari_column_vector( [complex_to_pari(z, working_prec) for z in manifold.tetrahedra_shapes('rect')] )


    manifold = manifold.copy()
    manifold.dehn_fill( (1, 0) ) 
    init_equations = manifold.gluing_equations('rect')

    

    target = sage_complex_to_pari(exp(target_meridian_log_holonomy), working_prec)
    
    if gluing_equation_error(init_equations, init_shapes, target) > pari(0.000001):
        raise ValueError('Initial solution not very good')

    # Now begin the actual computation
    
    eqns = enough_gluing_equations(manifold)
    assert eqns[-1] == manifold.gluing_equations('rect')[-1]
    
    shapes = init_shapes 
    for i in range(100):
        errors = gluing_equation_errors(eqns, shapes, target)
        if infinity_norm(errors) < target_espilon:
            break

        derivative = [ [  eqn[0][i]/z  - eqn[1][i]/(1 - z)  for i, z in enumerate(pari_vector_to_list(shapes))] for eqn in eqns]
        derivative[-1] = [ target*x for x in derivative[-1] ]
        derivative = pari_matrix(derivative)
        gauss = derivative.matsolve(pari_column_vector(errors))
        shapes = shapes - gauss

    # Check to make sure things worked out ok.
    error = gluing_equation_error(init_equations, shapes, target)
    total_change = infinity_norm(init_shapes - shapes)
    if error > 1000*target_espilon or total_change > pari(0.0000001):
        raise ValueError('Did not find a good solution to the gluing equations')


    manifold._cache["polished_shapes"] = shapes
    ans = pari_vector_to_list(shapes)
    if _within_sage:
        CC = ComplexField(bits_prec)
        ans = [CC(z) for z in ans]
    return ans


if __name__ == '__main__':
    import snappy
    M = snappy.Manifold('m071(7,0)')
    print polished_tetrahedra_shapes(M, (2*CC.pi()/7)*CC.gen(), bits_prec=10000)
    #M = snappy.Manifold('m071(0,0)')
    #print polished_tetrahedra_shapes(M,0, bits_prec=100)

    

    
