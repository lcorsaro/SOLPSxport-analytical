import numpy as np

from SOLPSxport import SOLPSxport
from SOLPSutils import read_dsa

from DChi import alphafunc_analytic,Darr_semianalytic,Chiarr_semianalytic,SOLPS2spline

def main(workdir,gfile_loc,lambdaq,n0,T0):
    SOLPSobj = SOLPSxport(workdir=workdir,gfile_loc=gfile_loc,impurity_list=[])

    SOLPSobj.getSOLPSlast10Profs()
    SOLPSobj.getSOLPSfluxProfs()
    dsa = SOLPSobj.data['solpsData']['dsa'] = read_dsa()
    GammaS=(SOLPSobj.data['solpsData']['profiles']['fluxD']-SOLPSobj.data['solpsData']['profiles']['fluxConv'])
    GammatotS=SOLPSobj.data['solpsData']['profiles']['fluxTot']
    qeS=SOLPSobj.data['solpsData']['profiles']['qe']
    qiS=SOLPSobj.data['solpsData']['profiles']['qi']
    TeS=SOLPSobj.data['solpsData']['last10']['te']
    TiS=SOLPSobj.data['solpsData']['last10']['ti']

    Lambdan=21/4
    LambdaT=7/2
    # nfunc=alphafunc_analytic(n0,lambdaq,Lambdan)
    # Tfunc=alphafunc_analytic(T0,lambdaq,LambdaT)
    # qfunc=alphafunc_analytic(q0,lambdaq)
    import pdb ; pdb.set_trace()
    lambdaq=float(lambdaq);n0=float(n0);T0=float(T0)
    dsa = np.array(dsa)
    D = Darr_semianalytic(lambdaq,n0,Lambdan,GammaS,dsa)
    Chii = Chiarr_semianalytic(lambdaq,n0,Lambdan,T0,LambdaT,GammatotS,qiS,TiS,dsa)
    Chie = Chiarr_semianalytic(lambdaq,n0,Lambdan,T0,LambdaT,GammatotS,qeS,TeS,dsa)

    SOLPSobj.data['solpsData']['xportCoef']={
        'dnew_flux':D,
        'kenew_flux':Chie,
        'kinew_flux':Chii,
        'vr_carbon':None, 'D_carbon':None,
    }

    write = False
    import pdb ; pdb.set_trace()

    if write:
        SOLPSobj.writeXport()

if __name__=='__main__':
    import sys
    _,workdir,gfile_loc,lambdaq=sys.argv
    main(workdir,gfile_loc,lambdaq,1e21,250)
    pass