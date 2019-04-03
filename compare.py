import numpy as np

mine = np.loadtxt("out_check")
ref  = np.loadtxt("../bead_modeling_flat/out_check")

mine_f = mine#[:,0]
ref_f = ref#[:,0]
err = []

check = ( len(mine_f) == len(ref_f) )
if( not check ):
    print("## LEN PROBLEM!! ##")
    print( f"REF: {len(ref_f)}, MINE: {len(mine_f)}" )
    exit()
else:
    print("Same len confirmed! ", len(ref_f), len(mine_f))

dim = min(mine.shape[0], ref.shape[0])
for i in range( dim ):
    rel_err = np.fabs( (mine_f[i] - ref_f[i])/ref_f[i] )
    err.append( rel_err )
    if( rel_err > 5e-4 ):
        print( f"CHECK! Rel err: {rel_err},  {mine_f[i]} != {ref_f[i]} at {i}" )

print( "Average rel err: ", np.mean(err) )
