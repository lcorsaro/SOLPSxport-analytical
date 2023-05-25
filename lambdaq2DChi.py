import numpy as np

from SOLPSxport import SOLPSxport
from SOLPSutils import read_dsa

from DChi import alphafunc_analytic,Darr_semianalytic,Chiarr_semianalytic,SOLPS2spline

def main(workdir,lambdaq):
    SOLPSobj = SOLPSxport(workdir=workdir,gfile_loc=None,impurity_list=[])

    SOLPSobj.getSOLPSlast10Profs()
    SOLPSobj.getSOLPSfluxProfs()
    dsa = SOLPSobj.data['solpsData']['dsa'] = read_dsa()
    GammaS=(SOLPSobj.data['solpsData']['profiles']['fluxD']-SOLPSobj.data['solpsData']['profiles']['fluxConv'])
    GammatotS=SOLPSobj.data['solpsData']['profiles']['fluxTot']
    qeS=SOLPSobj.data['solpsData']['profiles']['qe']
    qiS=SOLPSobj.data['solpsData']['profiles']['qi']
    TeS=SOLPSobj.data['solpsData']['last10']['te']
    TiS=SOLPSobj.data['solpsData']['last10']['ti']
    neS=SOLPSobj.data['solpsData']['last10']['ne']

    Lambdan=21/4
    LambdaT=7/2
    lambdaq=float(lambdaq)
    dsa = np.array(dsa)
    eV = 1.60217662e-19
    n0=SOLPS2spline(neS,dsa)(0)
    Te0=SOLPS2spline(TeS,dsa)(0)
    Ti0=SOLPS2spline(TiS,dsa)(0)
    D = Darr_semianalytic(lambdaq,n0,Lambdan,GammaS,dsa)
    Chii = Chiarr_semianalytic(lambdaq,n0,Lambdan,Ti0*eV,LambdaT,GammatotS,qiS,TiS*eV,dsa)
    Chie = Chiarr_semianalytic(lambdaq,n0,Lambdan,Te0*eV,LambdaT,GammatotS,qeS,TeS*eV,dsa)

    SOLPSobj.data['solpsData']['xportCoef']={
        'dnew_flux':D,
        'kenew_flux':Chie,
        'kinew_flux':Chii,
        'vr_carbon':None, 'D_carbon':None,
    }

    SOLPSobj.writeXport()

if __name__=='__main__':
    import sys
    _,workdir,lambdaq=sys.argv
    main(workdir,lambdaq)
    pass