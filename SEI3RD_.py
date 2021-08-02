#Autor: Juan Diego Mejía Becerra
#Correo: judmejiabe@unal.edu.co

#Libraries
import numpy as np
import plotly.graph_objs as go
import pandas as pd

from plotly.subplots import make_subplots
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

from scipy.integrate import odeint
from scipy.interpolate import interp1d


class SEI3RD(object):
    
    def __init__(self, 
                 incubationRate,
                 transmissionRate2,
                 transmissionRate3,
                 recoveryProbability0,
                 recoveryProbability1,
                 recoveryProbability2,
                 recoveryProbability3,
                 recoveryRate0,
                 recoveryRate1,
                 recoveryRate2,
                 recoveryRate3,
                 worsenCondition0,
                 worsenCondition1,
                 worsenCondition2,
                 worsenCondition3):
        """Crea una instancia de la clase SEIAR, que modelará una epidemia
        de acuerdo a un modelo SEI3RD con los parámetros dados"""
        self.omega = incubationRate
        self.beta2 = transmissionRate2
        self.beta3 = transmissionRate3
        self.delta0 = recoveryProbability0
        self.delta1 = recoveryProbability1
        self.delta2 = recoveryProbability2
        self.delta3 = recoveryProbability3
        self.gamma0 = recoveryRate0
        self.gamma1 = recoveryRate1
        self.gamma2 = recoveryRate2
        self.gamma3 = recoveryRate3
        self.sigma0 = worsenCondition0
        self.sigma1 = worsenCondition1
        self.sigma2 = worsenCondition2
        self.sigma3 = worsenCondition3
        
    def set_beta0(self, beta0):
        """Definir beta0(t)"""
        self.beta0 = beta0
        
    def set_beta1(self, beta1):
        """Definir beta1(t)"""
        self.beta1 = beta1
    
    def instantaneous_basic_reproductive_number(self, t):
        """Calcula el radio espectral de la matríz de la siguiente generación"""
        omega = self.omega
        beta0 = self.beta0
        beta1 = self.beta1
        beta2 = self.beta2
        beta3 = self.beta3
        delta0 = self.delta0
        delta1 = self.delta1
        delta2 = self.delta2
        delta3 = self.delta3
        gamma0 = self.gamma0
        gamma1 = self.gamma1
        gamma2 = self.gamma2
        gamma3 = self.gamma3
        sigma0 = self.sigma0
        sigma1 = self.sigma1
        sigma2 = self.sigma2
        sigma3 = self.sigma3
        
        self.T = np.matrix([[0, beta0(t), beta1(t), beta2, beta3],
                            [0, 0,     0,     0,     0],
                            [0, 0,     0,     0,     0],
                            [0, 0,     0,     0,     0],
                            [0, 0,     0,     0,     0]])
        
        #print(self.T)
        
        self.Sigma = np.matrix([[- omega, 0, 0, 0, 0],
                            [omega,   - (delta0 * gamma0 + (1 - delta0) * sigma0), 0, 0, 0],
                            [0, (1 - delta0) * sigma0, - (delta1 * gamma1 + (1 - delta1) * sigma1), 0, 0],
                            [0, 0, (1 - delta1) * sigma1, - (delta2 * gamma2 + (1 - delta2) * sigma2), 0],
                            [0, 0, 0, (1 - delta2) * sigma2, - (delta3 * gamma3 + (1 - delta3) * sigma3)]])
        
        return(np.max(np.linalg.eigvals(- np.matmul(self.T, np.linalg.inv(self.Sigma)))))
            
    
    def differential_equation(self, y, t):
        """Define el sistema de ecuaciones diferenciales a resolver"""
        omega = self.omega
        beta0 = self.beta0
        beta1 = self.beta1
        beta2 = self.beta2
        beta3 = self.beta3
        delta0 = self.delta0
        delta1 = self.delta1
        delta2 = self.delta2
        delta3 = self.delta3
        gamma0 = self.gamma0
        gamma1 = self.gamma1
        gamma2 = self.gamma2
        gamma3 = self.gamma3
        sigma0 = self.sigma0
        sigma1 = self.sigma1
        sigma2 = self.sigma2
        sigma3 = self.sigma3
        
        S, E, I0, I1, I2, I3, R, D = y
        dSdt = - S * (beta0(t) * I0 + beta1(t) * I1 + beta2 * I2 + beta3 * I3)
        dEdt = S * (beta0(t) * I0 + beta1(t) * I1 + beta2 * I2 + beta3 * I3) - omega * E
        dI0dt = omega * E - delta0 * gamma0 * I0 - (1 - delta0) * sigma0 * I0
        dI1dt = (1 - delta0) * sigma0 * I0 - delta1 * gamma1 * I1 - (1 - delta1) * sigma1 * I1
        dI2dt = (1 - delta1) * sigma1 * I1 - delta2 * gamma2 * I2 - (1 - delta2) * sigma2 * I2
        dI3dt = (1 - delta2) * sigma2 * I2 - delta3 * gamma3 * I3 - (1 - delta3) * sigma3 * I3
        dRdt = delta0 * gamma0 * I0 + delta1 * gamma1 * I1 + delta2 * gamma2 * I2 + delta3 * gamma3 * I3
        dDdt = (1 - delta3) * sigma3 * I3
        
        return(dSdt, dEdt, dI0dt, dI1dt, dI2dt, dI3dt, dRdt, dDdt)
    
    def solve_SEI3RD(self, E_0, I0_0, I1_0, I2_0, I3_0, R_0, D_0, T):
        """Resuelve las ecuaciones diferenciales y calcula el número final
        de individuos retirados. R_0 en este caso es la proporción inicial 
        de individuos retirados. T nos dá la ventana máxima de tiempo 
        donde se quiere observar el fenómeno"""
        self.E0 = E_0
        self.I0_0 = I0_0
        self.I1_0 = I1_0
        self.I2_0 = I2_0
        self.I3_0 = I3_0
        self.R_0 = R_0
        self.D_0 = D_0
        self.T = T
        
        self.t = np.linspace(0, T, 30000)
        
        S0 = 1 - E_0 - I0_0 - I1_0 - I2_0 - I3_0 - R_0 - D_0
        self.S0 = S0
        
        self.y0 = (S0, E_0, I0_0, I1_0, I2_0, I3_0, R_0, D_0)
        
        #Solución de las ecuaciones diferenciales
        ret = odeint(func = self.differential_equation, y0 = self.y0, t = self.t)
        
        self.S, self.E, self.I0, self.I1, self.I2, self.I3, self.R, self.D = ret.T
        
        self.S_ = interp1d(self.t, self.S)
        self.E_ = interp1d(self.t, self.E)
        self.I0_ = interp1d(self.t, self.I0)
        self.I1_ = interp1d(self.t, self.I1)
        self.I2_ = interp1d(self.t, self.I2)
        self.I3_ = interp1d(self.t, self.I3)
        self.R_ = interp1d(self.t, self.R)
        self.D_ = interp1d(self.t, self.D)
    
    def Rt_(self, t):
        R0t = []
        for t_ in self.t:
            R0t.append(self.instantaneous_basic_reproductive_number(t_))
        R0t = np.array(R0t)
        return(interp1d(self.t, R0t * self.S_(self.t))(t))
    
    def jc_table(self, numberSusceptible):
        
        seq = np.arange(0., 366.)

        results = pd.DataFrame({'Día' : seq,
                               'Susceptibles' : self.S_(seq) * numberSusceptible,
                               'Expuestos' : self.E_(seq) * numberSusceptible,
                               'Asintomáticos' : self.I0_(seq) * numberSusceptible,
                               'Moderados' : self.I1_(seq) * numberSusceptible,
                               'Severos' : self.I2_(seq) * numberSusceptible,
                               'Críticos' : self.I3_(seq) * numberSusceptible,
                               'Recuperados' : self.R_(seq) * numberSusceptible,
                               'Muertos' : self.D_(seq) * numberSusceptible,
                               'Rt' : self.Rt_(seq)})

        return(results)
        
    def plot_SEI3RD(self, numberSusceptible, separated = True):
        # Create traces
        
        w = 1.2
        
        trace0 = go.Scatter(
                x = self.t,
                y = self.S * numberSusceptible,
                mode = 'lines',
                name = 'Susceptibles',
                line = dict(color = 'blue',
                            width = w)
                )
        
        trace1 = go.Scatter(
                x = self.t,
                y = self.E * numberSusceptible,
                mode = 'lines',
                name = 'Expuestos',
                line = dict(color = 'pink',
                            width = w)
                )
        
        trace2 = go.Scatter(
                x = self.t,
                y = self.I0 * numberSusceptible,
                mode = 'lines',
                name = 'Asintomáticos',
                line = dict(color = 'cyan', 
                            width = w)
                )
        
        trace3 = go.Scatter(
                x = self.t,
                y = self.I1 * numberSusceptible,
                mode = 'lines',
                name = 'Síntomas Moderados',
                line = dict(color = 'yellow',
                            width = w)
                )
        
        trace4 = go.Scatter(
                x = self.t,
                y = self.I2 * numberSusceptible,
                mode = 'lines',
                name = 'Síntomas Severos',
                line = dict(color = 'orange',
                            width = w)
                )
        
        trace5 = go.Scatter(
                x = self.t,
                y = self.I3 * numberSusceptible,
                mode = 'lines',
                name = 'Críticos',
                line = dict(color = 'red',
                            width = w)
                )
        
        trace6 = go.Scatter(
                x = self.t,
                y = self.R * numberSusceptible,
                mode = 'lines',
                name = 'Recuperados',
                line = dict(color = 'green',
                            width = w)
                )
        
        trace7 = go.Scatter(
                x = self.t,
                y = self.D * numberSusceptible,
                mode = 'lines',
                name = 'Muertes',
                line = dict(color = 'black',
                            width = w)
                )
        
        if separated:
        
            fig = make_subplots(rows=4, cols=2,
            subplot_titles=("Susceptibles", 
                            "Expuestos (Estado Latente)", 
                            "Asintomáticos", 
                            "Síntomas Moderados",
                            "Síntomas Severos",
                            "Críticos",
                            "Recuperados",
                            "Muertos"))
        
            fig.add_trace(trace0, row = 1, col = 1)
            fig.add_trace(trace1, row = 1, col = 2)
            fig.add_trace(trace2, row = 2, col = 1)
            fig.add_trace(trace3, row = 2, col = 2)
            fig.add_trace(trace4, row = 3, col = 1)
            fig.add_trace(trace5, row = 3, col = 2)
            fig.add_trace(trace6, row = 4, col = 1)
            fig.add_trace(trace7, row = 4, col = 2)
            fig.update_layout(title = 'Proyecciones Modelo SEI3RD')
            
        else:
            fig = go.Figure()
            fig.add_trace(trace0)
            fig.add_trace(trace1)
            fig.add_trace(trace2)
            fig.add_trace(trace3)
            fig.add_trace(trace4)
            fig.add_trace(trace5)
            fig.add_trace(trace6)
            fig.add_trace(trace7)
            
            fig.update_layout(
                title="Proyecciones Modelo SEI3RD",
                xaxis_title = "Días",
                yaxis_title = "Prevalencia",
            )
        
        fig.layout.template = 'seaborn'
        fig.show()


        