from gassy.two_body_evolver import evolve_bodies
from scipy.io import loadmat
from deepdiff import DeepDiff


def test_evolver():
    """
    [X,Y,v_x,v_y,t] = spiralling113MESA15(100,10,1,0,0.5,10)
    out = struct('X', X, 'Y', Y, 'v_x', v_x, 'v_y', v_y, 't', t)
    save('spiral_out' out)
    """
    matlab_results = loadmat("tests/data/spiral_out.mat")

    print(matlab_results['out'])
    
    X,Y,v_x,v_y,t = evolve_bodies(M=100,m=10,a=1,e=0,dt=0.5,Tend=10, mesa_profile_name="15Msol")
    python_results = dict(X=X,Y=Y,v_x=v_x,v_y=v_y,t=t)
    diffs = DeepDiff(matlab_results, python_results)
    print(diffs)
    assert False
