import numpy as np
import pdb

class PreHeat:
    """
    FLAME SPEED AND STRAIN RATE AT PRE-HEAT ZONE
    Obtain the location of the max. strain rate upstream of the pre-heat zone
    and the extrapolated flame speed. 
    """
    def __init__(self,ref_Scalar,x,velocity,rho):
        """
        Reference plane:
        Any scalar to use as reference. Here, it was suggested the maximum value
        of O2 consumption
        """
        self.ref_Scalar = ref_Scalar
        self.refPoint = self.ref_Scalar.argmin()

        self.x = x
        self.velocity = velocity
        self.rho = rho

    def flameSpeed_Ku(self):
        Ku = np.gradient(self.velocity,self.x)
        maxStrLocation = abs(Ku).argmax()
        minVelocityPoint = self.velocity[:maxStrLocation].argmin()

        # S_mass = m(x)/rho_u
        S_mass_u = self.rho*self.velocity/self.rho[0]

        # Characteristic Strain Rate = K
        strainRatePoint = abs(Ku[:minVelocityPoint]).argmax()
        self.Ku_local = abs(Ku[strainRatePoint])

        # Characteristic Flame Speed Su
        Su_ = self.velocity[strainRatePoint] - self.Ku_local*(self.x[strainRatePoint:] - self.x[strainRatePoint])
        try:
            self.Su = Su_[self.refPoint]
        except:
            self.Su = np.nan
            self.Ku_local = np.nan
  
    def strain_rate_u(self):
        # characteristic Strain Rate = K
        self.flameSpeed_Ku()
        return self.Ku_local

    def flame_speed_u(self):
        """
        Su is taken as an extroplation of local strain rate
        """
        self.flameSpeed_Ku()
        return self.Su

class Reaction():
    """
    FLAME SPEED AND STRAIN RATE AT REACTION ZONE
    Obtained from the extrapolation of mass frow at burnt gas. 
    """

    def __init__(self,ref_Scalar,x,velocity,rho):
        """
        Reference plane:
        Any scalar to use as reference. Here, it was suggested the maximum value
        of O2 consumption
        """
        self.ref_Scalar = ref_Scalar
        self.refPoint = self.ref_Scalar.argmin()

        self.x = x
        self.velocity = velocity
        self.rho = rho

    def flameSpeed_Kb(self):

        # S_mass = m(x)/rho_b
        S_mass = self.rho*self.velocity/self.rho[self.rho.size-1]

        """ Kb = 1/rho*(grad(rho*u)) """
        grad = np.gradient(self.rho[self.refPoint:-1]*self.velocity[self.refPoint:-1],self.x[self.refPoint:-1])
        #Kb = np.array(-(1./(self.rho[self.refPoint:-1]*grad)))
        Kb = np.array(-(1./(self.rho[self.refPoint:-1]))*grad)
        self.Kb_local = Kb[0]

        # Weighted displacement speed
        self.Sl_d = self.rho[self.x.size-1]*S_mass[self.refPoint]/self.rho[0]

    def strain_rate_b(self):
        # characteristic Strain Rate = K
        self.flameSpeed_Kb()
        return self.Kb_local

    def flame_speed_b(self):
        """
        Sb is taken as an extroplation of local strain rate
        """
        self.flameSpeed_Kb()
        return self.Sl_d    

class FlameSpeeds:

    def __init__(self,x,rho,Y,wdot):
        self.x = x
        self.rho = rho
        self.Y = Y
        self.wdot = wdot
    
    def consumption_speed(self):

        Yb = self.Y[self.Y.size-1]
        Yu = self.Y[0]
        rho_u = self.rho[0]
        integral_wdot = np.trapz(self.wdot, self.x)
        Sc = abs(integral_wdot)/(Yu - Yb)/rho_u

        return Sc

