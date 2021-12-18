import matplotlib.pyplot as plt
from modsim import *
#pip install modsimpy and pint

def make_system(alpha, beta, gamma, delta, epsilon, zeta, eta, theta, iota, kappa, labda, mu, nu, xi, omicron, phi):
    """Make A system object for the SIR model.
    
    beta: contact rate in days
    gamma: recovery rate in days
    
    returns: System object
    """
    init = State(S=385600, I=2410, D=0, A=0, R=0, T=0, H=0, E=0) #s=20% or 385600 people susceptible, I=2410 people infected
    init /= sum(init) #convert from number ofpeople to fractions

    t0 = 0
    t_end = 7 * 36 #36 weeks

    return System(init=init, t0=t0, t_end=t_end,alpha = alpha,beta = beta,gamma = gamma,delta = delta,epsilon = epsilon,
                  theta = theta,eta = eta,zeta = zeta,mu = mu,nu = nu,tau = tau,lamda = lamda,kappa = kappa,xi = xi,rho = rho,sigma = sigma)

tc = 15      # time between contacts in days
tr = 14      # recovery time in days
#transmission rate
alpha = 0.4
beta = 0.3
gamma = 0.2
delta = 0.1
#rate of detection
epsilon = 0.5
theta = 0.6
#incubation rate
eta = 0.45
zeta = 0.88
#life threatening incubation rate
mu = 0.2
nu = 0.2
#mortality rate
tau = 0.03
#recovery rate
lamda = 0.6
kappa = 0.54
xi = 0.57
rho = 0.5
sigma = 0.4
# beta = 1 / tc      # contact rate in per day
# gamma = 1 / tr     # recovery rate in per day
system = make_system(alpha, beta, gamma, delta, epsilon, theta, eta, zeta, mu, nu, tau, lamda, kappa, xi, rho, sigma)

def update_func(state, T, system):
    """Update the SIR model.
    
    state: State with variables S, I, R
    T: time step
    system: System with beta and gamma
    
    returns: State object
    """
    constants  = ["alpha", "beta", "gamma", "delta", "epsilon", "theta", "eta", "zeta", "mu", "nu", "tau", "lamda", "kappa", "xi", "rho", "sigma"]
    (alpha, beta, gamma, delta, epsilon, theta, eta, zeta, mu, nu, tau, lamda, kappa, xi, rho, sigma) = [system[key] for key in constants]
    S,I,D,A,R,T,H,E = state

    infected = S*(alpha*I+beta*D+gamma*A*delta*R) - (epsilon + zeta + lamda )*I
    diagnosed = epsilon*I-(eta+rho)*D
    ailing = zeta*I-(theta+mu+kappa)*A
    recognized = eta*D+theta*A-(nu+tau)*R
    threatened = mu*A+nu*R-(sigma+tau)*T
    healed = lamda*I+rho*D+kappa*A+xi*R+sigma*T
    extinct = tau*T
    
    S -= infected
    I += infected
    D += diagnosed
    A += ailing 
    R += recognized 
    T += threatened 
    H += healed
    E += extinct
    return State(S=S,I=I,D=D,A=A,R=R,T=T,H=H,E=E)

init = State(S=385600, I=2410, D=0, A=0, R=0, T=0, H=0, E=0)
init /= sum(init)
state = update_func(init, 0, system)

def run_simulation(system, update_func):
    """Runs a simulation of the system.
    
    Add three Series objects to the System: S, I, R
    
    system: System object
    update_func: function that updates state
    """
    frame = TimeFrame(columns=system.init.index)
    frame.row[system.t0] = system.init

    state = system.init
    t0 = system.t0
    
    for t in linrange(system.t0, system.t_end):
        frame.row[t+1] = update_func(frame.row[t], t, system)
    
    return frame

tc = 9      # time between contacts in days, it can be 9 days already 
tr = 20     # recovery time in days
#transmission rate * contacts per day
alpha = .005*6
beta = .03*6
gamma = .003*9
delta = .005*9
#rate of detection
epsilon = 0.03/7
theta = 0.06*0.07
#incubation rate
eta = 0.2/7 
zeta = 0.2/6
#life threatening incubation rate
mu=0.4*.02
nu=.03/7
# mortality rate
tau = .03
#recovery rate
lamda = 0.4*0.01
kappa = 0.25*0.03*0.07
xi = 0.3*0.05*0.07
rho = 0.4*0.003
sigma = 0.05*0.003

system = make_system(alpha, beta, gamma, delta, epsilon, theta, eta, zeta, mu, nu, tau, lamda, kappa, xi, rho, sigma)
results = run_simulation(system, update_func)
results.head()

def plot_results(S, I, D, A, R, T, H,E):  
    plot(S, 'yellow', label='Susceptible')
    plot(I,'red', label = "Infected")
    plot(D,'magenta', label = "Diagnosed")
    plot(A,'purple', label = "Ailing")
    plot(R,'cyan', label = "Recognized")
    plot(T,'orange', label = "Threatened")
    plot(H,'green', label = "Healed")
    plot(E,'black', label = "Extinct")
    decorate(xlabel='Time (days)', ylabel='Population (1.0 = 100%)')                
plot_results(results.S, results.I, results.D, results.A, results.R, results.T, results.H,results.E)