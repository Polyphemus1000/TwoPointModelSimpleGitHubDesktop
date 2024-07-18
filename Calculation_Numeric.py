from ReducedCubicEquation_InitialGuess  import ReducedCubicEquation_InitialGuess
from ReducedCubicEquation_NewtonRaphson import ReducedCubicEquation_NewtonRaphson
import math

'''
This program solves a reduced cubic eqaation as set out in Fulvio Mitellos book page 268


'''

BoltzmannConstant= 1.38 * 10**(-23)
ElectricalCharge= 1.602 * 10**(-19)
NumberDensitySOLMid =  10**19
ParallelHeatFlux = 5 * 10**6
MassDeuterium = 3.3435 * 10**(-27)
Gamma = 7.5
TokamakMajorRadius = 9 # This is a very large tokamak
TokamakSafetyFactor = 3 # This is a large safety factor
ElectronParallelConductivityCoefficient=2000 # Common value that is used

'''
We use these default values which are higher than most tokamaks
'''

class Calculation_Numeric():


    def __init__(self):
        #print('Hello World')
        pass
    

    
    def Calculate(self, A, B, C, D, E):
    
        '''
        This function/method obtains the independent variable from a reduced cubic equation.
        The parameters for the reduced cubic equation are a, p and q i,e ax**3 + px + q 
        These are the values obtained by calculation in the first part of the method below
        '''
        #print('HelloAgain')
        try:
            NumberDensitySOLMid = A # not sure why you have to use *args for this to work
            ParallelHeatFlux = B # not sure why you have to use *args for this to work
            fpowerLoss = C
            fmomentumLoss = D
            fconductionLoss = E
            ConnectionLength =  math.pi * TokamakMajorRadius * TokamakSafetyFactor
            p = 3.5 * ParallelHeatFlux * ConnectionLength *fconductionLoss /ElectronParallelConductivityCoefficient
           # print(p)
            Beta = math.sqrt(2)* ElectricalCharge**(3.0/2)/math.sqrt(MassDeuterium)
           # print(Beta)
            qInitial = (2*ParallelHeatFlux*fpowerLoss/(Beta* NumberDensitySOLMid * fmomentumLoss))
            q= qInitial**(7.0/2)
            
            #print(q)
            
            a=1
            M= ReducedCubicEquation_InitialGuess(a,p,q)
            # The initial guess just involves coming up with a high and low value so we have bounds to apply the Newton Raphson
           # print('This is the reduced')
            #print(M.xlow)
            #print(M.xhigh)
            
            xlowhere= abs(M.xlow)
            xhighhere= abs(M.xhigh)
            
            P= ReducedCubicEquation_NewtonRaphson(xlowhere, xhighhere, a,p,q)
            self.TemperatureTarget = (P.xnew)**(4.0/7)
            
            
            
            TemperatureMidPontFirstPart = self.TemperatureTarget**(7.0/2) + p
            
            self.TemperatureMidPoint =  TemperatureMidPontFirstPart**(2.0/7)
            
            #print(self.TemperatureTarget) 
            #print(self.TemperatureMidPoint)
            self.NumberDensityTarget=  NumberDensitySOLMid * self.TemperatureMidPoint/(2 *  self.TemperatureTarget)
           # NumberDensityAtTarget.set(NumberDensityTarget)
           # print(self.NumberDensityTarget)
            self.SOLMidPointPressure = NumberDensitySOLMid * self.TemperatureMidPoint
            self.FluxAtTarget = self.NumberDensityTarget*math.sqrt(self.TemperatureTarget)* math.sqrt(2*ElectricalCharge/MassDeuterium)/self.SOLMidPointPressure
            
            self.TargetPressure = 2* self.NumberDensityTarget * self.TemperatureTarget
        except ValueError:
            pass
    
    
