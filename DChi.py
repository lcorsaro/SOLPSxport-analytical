import numpy as np
from scipy.interpolate import PchipInterpolator
from typing import Callable

def alphafunc_analytic(alpha0:float,lambdaq:float,Lambdaa:float=1)->Callable[[np.ndarray],np.ndarray]:
    def _alphafunc(rho:np.ndarray)->Callable[[np.ndarray],np.ndarray]:
        return alpha0*np.exp(-rho/(Lambdaa*lambdaq))
    return _alphafunc

def SOLPS2spline(arr:np.ndarray,rhoarr:np.ndarray):
    _SOLPS2spline = PchipInterpolator(x=rhoarr,y=arr,extrapolate=True)
    return _SOLPS2spline

def Dfunc_semianalytic(lambdaq:float,n0:float,Lambdan:float,
                       GammaS:Callable[[np.ndarray],np.ndarray])->Callable[[np.ndarray],np.ndarray]:
    n = alphafunc_analytic(n0,lambdaq,Lambdan)
    def _Dfunc(rho:np.ndarray)->Callable[[np.ndarray],np.ndarray]:
        return GammaS(rho)*Lambdan*lambdaq / n(rho)
    return _Dfunc

def Darr_semianalytic(lambdaq:float,n0:float,Lambdan:float,
                      GammaS:np.ndarray,
                      rho:np.ndarray)->np.ndarray:
    return (GammaS*Lambdan*lambdaq)*np.exp(rho/(lambdaq*Lambdan))/n0

def Chifunc_semianalytic(lambdaq:float,n0:float,Lambdan:float,T0:float,LambdaT:float,
                         GammaS:Callable[[np.ndarray],np.ndarray],
                         qS:Callable[[np.ndarray],np.ndarray],
                         TS:Callable[[np.ndarray],np.ndarray])->Callable[[np.ndarray],np.ndarray]:
    n = alphafunc_analytic(n0,lambdaq,Lambdan)
    T = alphafunc_analytic(T0,lambdaq,LambdaT)
    def _Chifunc(rho:np.ndarray)->Callable[[np.ndarray],np.ndarray]:
        return (qS(rho)-5/2 * GammaS(rho)*TS(rho))*LambdaT*lambdaq/(n(rho)*T(rho))
    return _Chifunc

def Chiarr_semianalytic(lambdaq:float,n0:float,Lambdan:float,T0:float,LambdaT:float,
                        GammaS:np.ndarray,qS:np.ndarray,TS:np.ndarray,
                        rho:np.ndarray)->np.ndarray:
    return LambdaT*lambdaq*(qS-5/2 * GammaS*TS)*np.exp(rho*(1/Lambdan+1/LambdaT)/lambdaq)/(n0*T0)

if __name__=='__main__':
    rho=np.linspace(0,5,100)
    Gamma=np.exp(-rho)
    n0=T0=q0=10
    lambdaq = 1.5
    Lambdan=21/4
    LambdaT=7/2
    nfunc=alphafunc_analytic(n0,lambdaq,Lambdan)
    Tfunc=alphafunc_analytic(T0,lambdaq,LambdaT)
    qfunc=alphafunc_analytic(q0,lambdaq)
    equalcheck = lambda arr1,arr2: np.all((np.abs(arr1-arr2)<=1e-15))
    Dfunc_=Dfunc_semianalytic(lambdaq,n0,Lambdan,SOLPS2spline(Gamma,rho))
    Darr_=Darr_semianalytic(lambdaq,n0,Lambdan,Gamma,rho)
    Chifunc_=Chifunc_semianalytic(lambdaq,n0,Lambdan,T0,LambdaT,SOLPS2spline(Gamma,rho),qfunc,Tfunc)
    Chiarr_=Chiarr_semianalytic(lambdaq,n0,Lambdan,T0,LambdaT,Gamma,qfunc(rho),Tfunc(rho),rho)
    print(equalcheck(Dfunc_(rho),Darr_))
    print(equalcheck(Chifunc_(rho),Chiarr_))